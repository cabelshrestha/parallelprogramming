/************************************************
CS789 - Parallel Programming 

	Implementation of FORWARDREDUCE, BACKWARDREDUCE
	and FASTRECUDE on a Hamiltonian cycle.

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
#define ORGSOURCE 6
#define REDUCE_STOP 22
#define VALUE 23

/* Message types */
#define MSG 0
#define ACK 1
#define DONE 2
#define STOP 3
#define PRINT 4
#define TREE 5
#define TREEDONE 6
#define PARENT 7
#define PARENTFOUND 8
#define FORBCAST 9
#define BCAST 10
#define FORBCASTACK 11
#define BCASTACK 12
#define REDUCEFORWARD 13
#define REDUCEBACKWARD 14
#define FASTREDUCE 15
#define INITDONE 16

/* Processes */
#define MASTER_PROCESS (size-1)
#define ROOT_NODE (size-2)
#define MEDIATOR 0

int D, M, rank, size, incoming, parent_found=0;
int i, j;
int next_dest, ack_dest, parent_dest;

int *msg, *ack_msg, *ack_forbcast;
bool reduceforward_done, reducebackward_done, fastreduce_done; 
bool reduce_initiated = false;
int send_count = 0;
int forbcast_count = 0;
int msg_count = 0;
int ack_count = 0;
int bcast_count = 0;
bool bcast_started = false;
bool bcast_ended = false;
bool init_done = false;

MPI_Status probe_status, status;


const char *type_map[] = {"MSG", "ACK", "DONE", "STOP", "PRINT", "TREE", "TREEDONE", "PARENT", "PARENTFOUND", "FORBCAST", "BCAST", "FORBCASTACK", "BCASTACK", "REDUCEFORWARD", "REDUCEBACKWARD", "FASTREDUCE", "INITDONE" };

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

unsigned int compute_parent_dest(unsigned int rank) {

	unsigned int parent = rank;
	unsigned int neighbor;

	for (i = 0; i < D; i++) {
		neighbor = rank ^ power2[i];
		if (neighbor > parent) {
			parent = neighbor;
		}
	}
	return parent;
}

void route_to_nextdestination(int* msg, int next_dest) {

	switch(msg[TYPE]) {
		case TREE:
		case STOP:
		case INITDONE:
			if (rank != MEDIATOR) { //Mediator is already added by the MASTER.
				msg[HOPCNT] = msg[HOPCNT] + 1;
				msg[HOPCNT + msg[HOPCNT] + 1] = rank;
			}	
			break;
		case PRINT:
		case DONE:
		case TREEDONE:
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

int msgCounter = 0;
void makeMessage(int *msg, int source, int dest, int type, int msgno) {
	msg[SOURCE] = source;
	msg[DEST] = dest;
	msg[TYPE] = type;
	msg[HOPCNT] = 0;
	// Done to retain original source for centralized printing.
	msg[HOPCNT+1] = source;
	if (msgno > 0)
		msg[MSGNO] = msgno;
	else {
		msg[MSGNO] = source*100+msgCounter;
		msgCounter++;
	}
}

int parity(int n) {
	int c = 0;
	while (n > 0) {
		if (n%2 == 1)
			c++;
		n >>= 1;
	}
	return c;
}

int findFirst(int i) {
	int c = 0;
	while (i % 2 == 0) {
		c++;
		i >>= 1;
	}
	return c;
}

int myNext(int i, int dim) {
	int c = 0; 

	/* For 000..00, we take first bit for swapping as it does not have any 1s. */
	if (i != 0) {
		c = findFirst(i);
	}

	int p = parity(i);
	int next;
	
	if (p%2 == 0) {
		next = i ^ power2[0];
	} else {
		if (i == power2[dim-1]) {
			next = 0;
		} else {
			next = i ^ power2[c+1];
		}
	}
	return next;
}

int myPrev(int i, int dim) {
	int c = 0; 

	/* For 000..00, we take first bit for swapping as it does not have any 1s. */
	if (i != 0) {
		c = findFirst(i);
	}

	int p = parity(i);
	int prev;
	
	if (p%2 == 0) {
		if (i == 0) {
			prev = power2[dim-1];
		} else {
			prev = i ^ power2[c+1];
		}
	} else {
		prev = i ^ power2[0];
		
	}
	return prev;
}

void initiate_reduce_forward(int i, int stop) {

	reduceforward_done = false;

	msg = (int *) malloc (sizeof (int) * MSGSIZE);

	makeMessage(msg, i, myNext(i, D), REDUCEFORWARD, 0);
	msg[ORGSOURCE] = msg[SOURCE];
	msg[VALUE] = rank+1;
	msg[REDUCE_STOP] = stop;
	MPI_Send(msg, MSGSIZE, MPI_INT, msg[DEST], 0, MPI_COMM_WORLD);
}

void initiate_reduce_backward(int i, int stop) {

	reducebackward_done = false;

	msg = (int *) malloc (sizeof (int) * MSGSIZE);
	makeMessage(msg, i, myPrev(i, D), REDUCEBACKWARD, 0);
	msg[ORGSOURCE] = msg[SOURCE];
	msg[VALUE] = rank+1;
	msg[REDUCE_STOP] = stop;
	MPI_Send(msg, MSGSIZE, MPI_INT, msg[DEST], 0, MPI_COMM_WORLD);
		
}

void printSummaryAtDestination(int *msg) {
	char buffer1[512], buffer2[128];

	sprintf(buffer1, "%d(%s): Message (%s) #%03d from %d(%s) to %d(%s) in %d hops ", msg[SOURCE], byte2bin(msg[SOURCE], D), type_map[msg[ORIGTYPE]], msg[MSGNO], msg[HOPCNT+1], byte2bin(msg[HOPCNT+1], D), msg[SOURCE], byte2bin(msg[SOURCE], D), msg[HOPCNT]);

	sprintf(buffer2, "(%s->", byte2bin(msg[HOPCNT+1], D));
	strcat(buffer1, buffer2);
	for (i = 0, j=1; i<msg[HOPCNT]; i++, j++) {
		sprintf(buffer2, "%s->", byte2bin(msg[HOPCNT+j+1], D));
		strcat(buffer1, buffer2);
	}
	sprintf(buffer2, "%s)", byte2bin(msg[SOURCE], D));
	strcat(buffer1, buffer2);

	if (msg[ORIGTYPE] == REDUCEFORWARD || msg[ORIGTYPE] == REDUCEBACKWARD) {
		sprintf(buffer2, " Reduce result=%d", msg[VALUE]);
		strcat(buffer1, buffer2);
	}

	if (msg[ORIGTYPE] == FASTREDUCE) {
		sprintf(buffer2, " Fastreduce result=%d", msg[VALUE]);
		strcat(buffer1, buffer2);
	}

	printf("%s\n", buffer1);
}


bool isMsgDone() {
	if (ack_count == msg_count)
		return true;
	return false;
}

bool isBroadcastDone() {
	if (bcast_ended)
		return true;
	return false;
}

bool isReduceForwardDone() {
	if (reduceforward_done)
		return true;
	return false;
}

bool isReduceBackwardDone() {
	if(reducebackward_done)
		return true;
	return false;
}

bool isFastReduceDone() {
	if (fastreduce_done)
		return true;
	return false;
}

bool isDone() {
	if ( 
		isMsgDone() 
		&& isBroadcastDone()
		&& isReduceForwardDone()
		&& isReduceBackwardDone()
		&& isFastReduceDone()
		) {
		return true;
	}

	return false;
}

void send_done_message() {
	msg = (int *) malloc (sizeof (int) * MSGSIZE);
	makeMessage(msg, rank, MEDIATOR, DONE, 0);
	next_dest = compute_next_dest(rank, msg[DEST]);
	MPI_Send(msg, MSGSIZE, MPI_INT, next_dest, 0, MPI_COMM_WORLD);
}


void master_code() {


	msg = (int *) malloc (sizeof (int) * MSGSIZE);

	/* Send TREE msg to all nodes to invoke spanning tree construction. */
	for(i = 0; i < size-1; i++) {
		printf("M: Injecting TREE msg for node-%d through mediator(node-0).\n", i);
		makeMessage(msg, MEDIATOR, i, TREE, 0);
		MPI_Send(msg, MSGSIZE, MPI_INT, MEDIATOR, 0, MPI_COMM_WORLD);
	}

	int done = 0;
	int treedone = false;

	while (done != pow(2, D) || !treedone) {

		MPI_Recv(msg, MSGSIZE, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG,MPI_COMM_WORLD, &status);

		if (msg[TYPE] == PRINT) {
			 printSummaryAtDestination(msg);
		} else if (msg[TYPE] == DONE) {
			done++;
		} else if (msg[TYPE] == TREEDONE) {
			printf("M: Spanning tree construction complete.\n");
			treedone = true;

			for(i = 0; i < size-1; i++) {
				printf("M: Injecting INITDONE msg for node-%d through mediator(node-0).\n", i);
				makeMessage(msg, MEDIATOR, i, INITDONE, 0);
				MPI_Send(msg, MSGSIZE, MPI_INT, MEDIATOR, 0, MPI_COMM_WORLD);
			}
		} 
	}

	printf("M: Injecting STOP msg using REDUCEFORWARD through mediator(node-0).\n");
	initiate_reduce_forward(MEDIATOR, 1);
}

void slave_code() {
	int children[D];
	int child_count = 0;
	int v = rank * 2;
	int reduce_count = 0;

	fastreduce_done = false;
	reduceforward_done = true;
	reducebackward_done = true;

	msg = (int *) malloc (sizeof (int) * MSGSIZE);
	int stop = false;

	if (rank%2 == 0) {
		if (rank == 0) {
			reduce_count = D;
		} else {
			reduce_count = findFirst(rank);
		}
	}

	while(!stop) {

		/* Probe to see if there is any incoming message. */
		MPI_Iprobe(MPI_ANY_SOURCE,0,MPI_COMM_WORLD,&incoming,&probe_status);

		if(incoming) {

			MPI_Recv(msg, MSGSIZE, MPI_INT, probe_status.MPI_SOURCE, 0, MPI_COMM_WORLD, &status);

			next_dest = compute_next_dest(rank, msg[DEST]);

			if (rank != next_dest && (msg[TYPE] != REDUCEFORWARD && msg[TYPE] != REDUCEBACKWARD)) { 
				/* If msg not for me, route it. */
				route_to_nextdestination(msg, next_dest);

			} else if ( rank != msg[ORGSOURCE] && (msg[TYPE] == REDUCEFORWARD || msg[TYPE] == REDUCEBACKWARD)) { //REDUCE* routing logic
				switch(msg[TYPE]) {
					case REDUCEFORWARD:
						msg[HOPCNT] = msg[HOPCNT] + 1;
						msg[HOPCNT + msg[HOPCNT] + 1] = rank;
						msg[SOURCE] = rank;
						msg[DEST] = myNext(rank, D);
						msg[VALUE] += rank+1;

						MPI_Send(msg, MSGSIZE, MPI_INT, msg[DEST], 0, MPI_COMM_WORLD);

						/* Terminate logic using reduce forward.*/
						if (msg[REDUCE_STOP] == 1)
							stop = true;

						break;

					case REDUCEBACKWARD:
						msg[HOPCNT] = msg[HOPCNT] + 1;
						msg[HOPCNT + msg[HOPCNT] + 1] = rank;
						msg[SOURCE] = rank;
						msg[DEST] = myPrev(rank, D);
						msg[VALUE] *= rank+1;

						MPI_Send(msg, MSGSIZE, MPI_INT, msg[DEST], 0, MPI_COMM_WORLD);

						/* Terminate logic using reduce backward.*/
						if (msg[REDUCE_STOP] == 1)
							stop = true;

						break;

				}

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

						if(isDone()) {
							send_done_message();
						}

						break;

					case TREE:

						/* Print msg. */
						send_for_printing(msg);

						parent_dest = compute_parent_dest(rank);

						/* Send parent msg. */
						if (parent_dest != rank) {
							/* This is only for non ROOT nodes.*/
							msg = (int *) malloc (sizeof (int) * MSGSIZE);
							makeMessage(msg, rank, parent_dest, PARENT, 0);
							next_dest = compute_next_dest(rank, msg[DEST]);
							MPI_Send(msg, MSGSIZE, MPI_INT, next_dest, 0, MPI_COMM_WORLD);
						}

						break;

					case PARENT:

						/* Log child. */
						children[child_count] = msg[SOURCE];
						child_count++;

						/* Print msg. */
						send_for_printing(msg);

						/* The ROOT does not need to send anything to itself.*/
						if ( rank != ROOT_NODE) {
							msg = (int *) malloc (sizeof (int) * MSGSIZE);
							makeMessage(msg, rank, ROOT_NODE, PARENTFOUND, 0);
							next_dest = compute_next_dest(rank, msg[DEST]);
							MPI_Send(msg, MSGSIZE, MPI_INT, next_dest, 0, MPI_COMM_WORLD);
						}	

						break;

					case PARENTFOUND: /* This case is only for ROOT node. */

						/* Print msg. */
						send_for_printing(msg);

						parent_found++;

						/* When tree is complete, send treedone to master through mediator.*/
						if (parent_found == pow(2, D) - (D+1)) {
							msg = (int *) malloc (sizeof (int) * MSGSIZE);
							makeMessage(msg, rank, MEDIATOR, TREEDONE, 0);
							next_dest = compute_next_dest(rank, msg[DEST]);
							MPI_Send(msg, MSGSIZE, MPI_INT, next_dest, 0, MPI_COMM_WORLD);
						}

						break;

					case PRINT:
					case DONE:
					case TREEDONE:

						/* This case is only for the MEDIATOR. Send everything to master. */
						MPI_Send(msg, MSGSIZE, MPI_INT, MASTER_PROCESS, 0, MPI_COMM_WORLD);

						break;

					case INITDONE:
						init_done = true;

						/* Print msg. */
						send_for_printing(msg);

						printf("%d(%s): Children are {", rank, byte2bin(rank, D));
						for (i = 0; i < child_count; i++) {
							printf("%d", children[i]);
							if (i != child_count-1)
								printf(", ");
						}
						printf("}\n");

						break;

					case FORBCAST: //ROOT_NODE receives this.
						
						/* Print msg. */
						send_for_printing(msg);
						
						bcast_started=true;
						/* Broadcast the message. */
						for(i = 0; i < child_count; i++) {
							msg = (int *) malloc (sizeof (int) * MSGSIZE);
							makeMessage(msg, ROOT_NODE, children[i], BCAST, 0);
							MPI_Send(msg, MSGSIZE, MPI_INT, msg[DEST], 0, MPI_COMM_WORLD); 
							bcast_count++;
						}

						break;

					case BCAST: //Nodes other than ROOT_NODE receives this.

						/* Print msg. */
						send_for_printing(msg);

						bcast_started = true;

						if (child_count > 0 ) {
							for(i = 0; i < child_count; i++) {
								msg = (int *) malloc (sizeof (int) * MSGSIZE);
								makeMessage(msg, rank, children[i], BCAST, 0);
								MPI_Send(msg, MSGSIZE, MPI_INT, msg[DEST], 0, MPI_COMM_WORLD); 
								bcast_count++;
							}
						} else {
							/* Send acknowledgement. */
							parent_dest = compute_parent_dest(rank);
							msg = (int *) malloc (sizeof (int) * MSGSIZE);
							makeMessage(msg, rank, parent_dest, BCASTACK, 0);
							MPI_Send(msg, MSGSIZE, MPI_INT, parent_dest, 0, MPI_COMM_WORLD);

							if(bcast_started)
								bcast_ended = true;

							if(isDone()) {
								send_done_message();
							}	
						}
						
						break;

					case FORBCASTACK:
						send_for_printing(msg);

						forbcast_count--;

						if(isDone()) {
							send_done_message();
						}	

						break;

					case BCASTACK:
						send_for_printing(msg);

						bcast_count--;

						if (bcast_count == 0) {
							if (rank == ROOT_NODE) {
								/* Send acknowledgement. */
								ack_forbcast = (int *) malloc (sizeof (int) * MSGSIZE);
								makeMessage(ack_forbcast, rank, MEDIATOR, FORBCASTACK, 0);
								next_dest = compute_next_dest(rank, ack_forbcast[DEST]);
								MPI_Send(ack_forbcast, MSGSIZE, MPI_INT, next_dest, 0, MPI_COMM_WORLD);	
							} else {
								/* Send acknowledgement. */
								parent_dest = compute_parent_dest(rank);

								msg = (int *) malloc (sizeof (int) * MSGSIZE);
								makeMessage(msg, rank, parent_dest, BCASTACK, 0);
								MPI_Send(msg, MSGSIZE, MPI_INT, parent_dest, 0, MPI_COMM_WORLD);
							}

							if(bcast_started)
								bcast_ended = true;

							if(isDone()) {
								send_done_message();
							}
						}

						break;

					case REDUCEFORWARD:

						reduceforward_done = true;

						if (msg[REDUCE_STOP] == 1) {
							stop = true;
						} else {
							send_for_printing(msg);

							if (ack_count == M && forbcast_count == 0 && bcast_count == 0 && reduceforward_done && reducebackward_done) {
								msg = (int *) malloc (sizeof (int) * MSGSIZE);
								makeMessage(msg, rank, MEDIATOR, DONE, 0);

								next_dest = compute_next_dest(rank, msg[DEST]);

								MPI_Send(msg, MSGSIZE, MPI_INT, next_dest, 0, MPI_COMM_WORLD);
							}
						}

						break;

					case REDUCEBACKWARD:
						reducebackward_done = true;

						if (msg[REDUCE_STOP] == 1) {
							stop = true;
						} else {
							send_for_printing(msg);

							if (ack_count == M && forbcast_count == 0 && bcast_count == 0 && reduceforward_done && reducebackward_done) {
								msg = (int *) malloc (sizeof (int) * MSGSIZE);
								makeMessage(msg, rank, MEDIATOR, DONE, 0);

								next_dest = compute_next_dest(rank, msg[DEST]);

								MPI_Send(msg, MSGSIZE, MPI_INT, next_dest, 0, MPI_COMM_WORLD);
							}

						}

						break;

					case FASTREDUCE:
						v += msg[VALUE];

						reduce_count--;

						if (reduce_count == 0) {
							if (rank != 0) {
								msg = (int *) malloc (sizeof (int) * MSGSIZE);
								makeMessage(msg, rank, rank ^ power2[findFirst(rank)], FASTREDUCE, 0);
								/* Setting a value for addition purpose. */
								msg[VALUE] = v;

								next_dest = compute_next_dest(rank, msg[DEST]);
								MPI_Send(msg, MSGSIZE, MPI_INT, next_dest, 0, MPI_COMM_WORLD);

							} else { //Fastreduce is done at this point.
								msg[VALUE] = v;
							}
							fastreduce_done = true;

							send_for_printing(msg);

							if(isDone()) {
								send_done_message();
							}
							 
						}

						break;

					case STOP:  //This is not being used for this question. Reduceforward is used for termination.

						stop = true;

						break;

				}
			}

		} else {
			if (send_count < M && init_done) {

				msg = (int *) malloc (sizeof (int) * MSGSIZE);

				/* 000...00 broadcasts its last msg instead of sending it to someone specific.*/
				if (rank == 0 && send_count == M-1) { //Send last msg of 0 for broadcast.
					makeMessage(msg, rank, ROOT_NODE, FORBCAST, 0);
					next_dest = compute_next_dest(rank, msg[DEST]);
					MPI_Send(msg, MSGSIZE, MPI_INT, next_dest, 0, MPI_COMM_WORLD); 
					forbcast_count++;
				} else {
					makeMessage(msg, rank, rand()%(size-1), MSG, 0);
					next_dest = compute_next_dest(rank, msg[DEST]);
					MPI_Send(msg, MSGSIZE, MPI_INT, next_dest, 0, MPI_COMM_WORLD); 
					msg_count++;		
				}

				send_count++;		

			}

			/* To test REDUCE* initiate forward and backward once from node-0. */
			if (!reduce_initiated && rank==0) {
				initiate_reduce_forward(rank, 0);
				initiate_reduce_backward(rank, 0);
				reduce_initiated = true;
			}

			/* Testing fast reduce implementation. */
			if ( (rank % 2)==1 && !fastreduce_done) {
				fastreduce_done = true;

				msg = (int *) malloc (sizeof (int) * MSGSIZE);
				makeMessage(msg, rank, rank^power2[0], FASTREDUCE, 0);
				/* Setting a value for addition purpose. */
				msg[VALUE] = v;

				next_dest = compute_next_dest(rank, msg[DEST]);
				MPI_Send(msg, MSGSIZE, MPI_INT, next_dest, 0, MPI_COMM_WORLD); 
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

	if (rank == MASTER_PROCESS) { 
		master_code();
	} else {
		slave_code();
	}

	MPI_Finalize();

	return 0;

}






