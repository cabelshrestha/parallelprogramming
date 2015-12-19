/************************************************
CS789 - Parallel Programming 

	Hypercube creation and routing.

@Author Cabel Dhoj Shrestha
Dec 2015
**************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <mpi.h>
#include <math.h>
#include <assert.h>
#include <sys/time.h>

#define assert__(x) for ( ; !(x) ; assert(x) )

#define MSGSIZE 25

/* Message index */
#define SOURCE 0
#define DEST 1
#define MSGNO 2
#define TYPE 3
#define ORIGTYPE 4
#define HOPCNT 5

/* Message types */
#define MSG 0
#define ACK 1
#define DONE 2
#define STOP 3
#define PRINT 4

/* Processes */
#define MASTER (size-1)
#define MEDIATOR 0

int D, M, rank, size, incoming;
int i, j;
int next_dest, ack_dest;
MPI_Status probe_status, status;
int *msg, *ack_msg;

const char *type_map[] = {"MSG", "ACK", "DONE", "STOP", "PRINT" };

int msgCounter = 0;
void makeMessage(int *msg, int source, int dest, int type, int msgno) {
	msg[SOURCE] = source;
	msg[DEST] = dest;
	msg[TYPE] = type;
	msg[HOPCNT] = 0;
	/* Done to retain original source for centralized printing.*/
	msg[HOPCNT+1] = source;
	if (msgno > 0)
		msg[MSGNO] = msgno;
	else {
		msg[MSGNO] = source*100+msgCounter;
		msgCounter++;
	}
}

const unsigned int power2[12] = {1,2,4,8,16,32,64,128,256,512,1024,2048};

unsigned int compute_next_dest(unsigned int rank, unsigned int dest) {
	if (rank == dest) // if we are already there
		return dest;
	unsigned int m = rank ^ dest;
	int c = 0;
	// find the first 1 searching from right to left.
	while (m%2 == 0) { // while m is even.
		c++;
		m >>= 1;
	}
	return power2[c] ^ rank;
}

void route_to_nextdestination(int* msg, int next_dest) {

	switch(msg[TYPE]) {
		case STOP:
			if (rank != MEDIATOR) { //Mediator is already added by the MASTER.
				msg[HOPCNT] = msg[HOPCNT] + 1;
				msg[HOPCNT + msg[HOPCNT] + 1] = rank;
			}	
			break;
		case PRINT:
		case DONE:
			break;
		default:
			msg[HOPCNT] = msg[HOPCNT] + 1;
			msg[HOPCNT + msg[HOPCNT] + 1] = rank;

	}

	MPI_Send(msg, MSGSIZE, MPI_INT, next_dest, 0, MPI_COMM_WORLD);

}

void send_for_printing(int* msg) {
	/* Send msg for centralized printing. */
	msg[ORIGTYPE] = msg[TYPE];
	msg[TYPE] = PRINT;
	msg[SOURCE] = rank;
	msg[DEST] = MEDIATOR;
	next_dest = compute_next_dest(rank, msg[DEST]);
	MPI_Send(msg, MSGSIZE, MPI_INT, next_dest, 0, MPI_COMM_WORLD);
} 

char *byte2bin(unsigned int x, int dim) {
	char* b = (char*) malloc(sizeof(char)*(dim));
	b[0] = '\0';
	int z;
	for (z = power2[dim-1]; z > 0; z >>= 1) {
		strcat(b, ((x & z) == z) ? "1" : "0");
	}
	return b;
}

void printSummaryAtDestination(int *msg) {
	char buffer1[128], buffer2[128];

	sprintf(buffer1, "%d(%s): Message (%s) #%03d from %d(%s) to %d(%s) in %d hops ", msg[SOURCE], byte2bin(msg[SOURCE], D), type_map[msg[ORIGTYPE]], msg[MSGNO], msg[HOPCNT+1], byte2bin(msg[HOPCNT+1], D), msg[SOURCE], byte2bin(msg[SOURCE], D), msg[HOPCNT]);
	sprintf(buffer2, "(%s->", byte2bin(msg[HOPCNT+1], D));
	strcat(buffer1, buffer2);
	for (i = 0, j=1; i<msg[HOPCNT]; i++, j++) {
		sprintf(buffer2, "%s->", byte2bin(msg[HOPCNT+j+1], D));
		strcat(buffer1, buffer2);
	}
	sprintf(buffer2, "%s)", byte2bin(msg[SOURCE], D));
	strcat(buffer1, buffer2);
	printf("%s\n", buffer1);
}


void master_code() {

	msg = (int *) malloc (sizeof (int) * MSGSIZE);
	int done = 0;

	// Receive a DONE and PRINT message from each worker 
	while (done != (size-1)) {

		MPI_Recv(msg, MSGSIZE, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG,MPI_COMM_WORLD, &status);

		if (msg[TYPE] == PRINT) {
			 printSummaryAtDestination(msg);
		} else if (msg[TYPE] == DONE) {
			done++;	
		}
	}

	// Send a STOP message to the MEDIATER (rank=0) for it to pass it on to other nodes.
	for (i =size-2; i>=0; i--) {
		makeMessage(msg, MEDIATOR, i, STOP, 0);
		MPI_Send(msg, MSGSIZE, MPI_INT, MEDIATOR, 0, MPI_COMM_WORLD);
	}
}

void slave_code() {
	int send_count = 0;
	int ack_count = 0;
	bool stop = false;

	while(!stop) {
		/* Probe to see if there is any incoming message. */
		MPI_Iprobe(MPI_ANY_SOURCE,0,MPI_COMM_WORLD,&incoming,&probe_status);

		if(incoming) {

			MPI_Recv(msg, MSGSIZE, MPI_INT, probe_status.MPI_SOURCE, 0, MPI_COMM_WORLD, &status);

			next_dest = compute_next_dest(rank, msg[DEST]);

			if (rank != next_dest) { //If msg not for me, route it.
				route_to_nextdestination(msg, next_dest);
			} else {
				switch(msg[TYPE]) {
					case MSG:
						/* Send the acknowledgement. */
						ack_dest = msg[SOURCE];	
						ack_msg = (int *) malloc (sizeof (int) * MSGSIZE);
						makeMessage(ack_msg, rank, ack_dest, ACK, 0);
						next_dest = compute_next_dest(rank, ack_msg[DEST]);
						MPI_Send(ack_msg, MSGSIZE, MPI_INT, next_dest, 0, MPI_COMM_WORLD);	

						/* Print msg. */
						send_for_printing(msg);
						
						break;

					case ACK:

						/* Print msg. */
						send_for_printing(msg);

						ack_count++;

						if (ack_count == M) {
							msg = (int *) malloc (sizeof (int) * MSGSIZE);
							makeMessage(msg, rank, MEDIATOR, DONE, 0);

							next_dest = compute_next_dest(rank, msg[DEST]);

							MPI_Send(msg, MSGSIZE, MPI_INT, next_dest, 0, MPI_COMM_WORLD);
						}
						break;
						
					case PRINT:
					case DONE:
						MPI_Send(msg, MSGSIZE, MPI_INT, MASTER, 0, MPI_COMM_WORLD);
						break;

					case STOP:
						if (next_dest == rank) {
							stop = true;
							break;
						} else {
							MPI_Send(msg, MSGSIZE, MPI_INT, next_dest, 0, MPI_COMM_WORLD); 
						}
						break;
				}

			}

		} else {
			if (send_count < M) {
				msg = (int *) malloc (sizeof (int) * MSGSIZE);
				makeMessage(msg, rank, rand()%(size-1), MSG, 0);

				next_dest = compute_next_dest(rank, msg[DEST]);

				MPI_Send(msg, MSGSIZE, MPI_INT, next_dest, 0, MPI_COMM_WORLD); 
				send_count++;	
			}
		}
	}

}

int main(int argc, char *argv[ ]) {

	if (argc != 3) {
		printf("Usage:\n hypercube <dimension> <# of messages>\n");
		exit(1);
	}

	D = atoi(argv[1]);
	M = atoi(argv[2]);

	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	MPI_Comm_size (MPI_COMM_WORLD, &size);

	/* Check necessary number of processes.*/
	assert__((size) == pow(2, D)+1) {
		if (rank == 1)
			printf("Error: Number of processes should be 2^d+1=2^%d+1=%.0f.\n", D, pow(2, D)+1);
		MPI_Finalize();
		exit(0);

	}

	struct timeval t;
	gettimeofday(&t,NULL);
	srand(t.tv_usec);

	if (rank == MASTER) { // Master
		master_code();
	} else {
		slave_code();
	}

	MPI_Finalize();

	return 0;

}






