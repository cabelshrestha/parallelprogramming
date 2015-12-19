/************************************************
CS789 - Parallel Programming 

  Optimized pipeline implementation.

@Author: Cabel Dhoj Shrestha
Sept 2015
**************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <unistd.h>
#include <sys/time.h>
#include <stdarg.h>

#define TAG 0
#define END_TAG 99

int rank, size, incoming;
long  val;
MPI_Status status;
MPI_Request req;
struct timeval t0, t1, t2, t3;

void proc(int delay, int from, int to) {
  while(1) {
    MPI_Recv(&val, 1, MPI_LONG, from, TAG, MPI_COMM_WORLD, &status);
    usleep(delay*10000);
    MPI_Ssend(&val, 1, MPI_LONG, to,   TAG, MPI_COMM_WORLD);

    MPI_Iprobe(MPI_ANY_SOURCE,END_TAG,MPI_COMM_WORLD,&incoming,&status);
    if (incoming) {
      printf("%d is terminating!\n", rank);
      MPI_Recv(&val, 1, MPI_INT, status.MPI_SOURCE, END_TAG, MPI_COMM_WORLD, &status);
      break;
    }
  }
}

void in(int no, int to) {
  int i;
  for(i=0; i<no; i++) {
    gettimeofday(&t1,NULL);
    val = (t1.tv_sec%1000)*1000000+t1.tv_usec;    
    MPI_Ssend(&val, 1, MPI_LONG, to, TAG, MPI_COMM_WORLD);
  }
}
  
void out(int no, int from) {
  int i;
  long long elapsedsum =0;
  long elapsed;
  long long transportsum = 0;
  long transport;
  int getlatency = 1;
  long latency;
  gettimeofday(&t0,NULL);
  gettimeofday(&t3,NULL);
  for(i=0; i<no; i++) {
    MPI_Recv(&val, 1, MPI_LONG, from, TAG, MPI_COMM_WORLD, &status);
    gettimeofday(&t2,NULL);

    elapsed = (t2.tv_sec%1000)*1000000+t2.tv_usec-(t3.tv_sec%1000)*1000000-t3.tv_usec;
    if (getlatency) {
      latency = elapsed;
      getlatency = 0;
    }
    else 
      elapsedsum += elapsed;
    transport = (t2.tv_sec%1000)*1000000+t2.tv_usec-val;
    transportsum += transport;
    printf("Packet %d - Transport Time ..........: %ld   ", i, transport);
    printf("Elapsed since last packet: %ld\n", elapsed);    
    gettimeofday(&t3,NULL);
  }
  gettimeofday(&t2,NULL);
  
  printf("Latency...............: %ld\n",latency);
  printf("Total Time............: %ld\n",(t2.tv_sec%1000)*1000000+t2.tv_usec-(t0.tv_sec%1000)*1000000-t0.tv_usec);

  printf("Average Transport Time: %ld\n", (long)(transportsum/no));
  printf("Average Rate..........: %ld\n", (long)(elapsedsum/(no-1)));

  printf("Press enter to terminate\n");
  getchar();

  /* Send termination message. */
  for (i=0;i<size;i++) {
      MPI_Send(&val, 1, MPI_INT, i, END_TAG, MPI_COMM_WORLD);
  }
}


void disperse(int size, int from, ...) {
  va_list argp;
  int i, incoming;
  while(1) {
    va_start(argp, from);
    for (i=0;i<size;i++) {
      MPI_Recv(&val, 1, MPI_LONG, from, TAG, MPI_COMM_WORLD, &status);
      MPI_Ssend(&val, 1, MPI_LONG, va_arg(argp, int), TAG, MPI_COMM_WORLD);  
    }
    va_end(argp);

    MPI_Iprobe(MPI_ANY_SOURCE,END_TAG,MPI_COMM_WORLD,&incoming,&status);
    if (incoming) {
      printf("%d is terminating!\n", rank);
      MPI_Recv(&val, 1, MPI_INT, status.MPI_SOURCE, END_TAG, MPI_COMM_WORLD, &status);
      break;
    }

  }
}


void collect(int size, int to, ...) {
  va_list argp;
  int i, incoming=0;
  while(1) {
    va_start(argp, to);
    for (i=0;i<size;i++) {
      MPI_Recv(&val, 1, MPI_LONG, va_arg(argp, int), TAG, MPI_COMM_WORLD, &status);
      MPI_Ssend(&val, 1, MPI_LONG, to, TAG, MPI_COMM_WORLD);
    }
    va_end(argp);

    MPI_Iprobe(MPI_ANY_SOURCE,END_TAG,MPI_COMM_WORLD,&incoming,&status);
    if (incoming) {
      printf("%d is terminating!\n", rank);
      MPI_Recv(&val, 1, MPI_INT, status.MPI_SOURCE, END_TAG, MPI_COMM_WORLD, &status);
      break;
    }
  }
}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int no = atoi(argv[1]);

  switch(rank) {
    case 0: in(no, 1);  break;
    case 1: proc(2, 0, 2); break;
    case 2: disperse(3, 1, 3,4,5); break;
    case 3:  
    case 4:  
    case 5: proc(6, 2, 6); break;
    case 6: collect(3, 7, 3,4,5); break;
    case 7: disperse (3, 6, 8, 9, 10); break;
    case 8: 
    case 9: 
    case 10: proc(6, 7, 11); break;
    case 11: collect(3, 12, 8,9,10); break;
    case 12: disperse (4, 11, 13,14,15,16); break;
    case 13: 
    case 14: 
    case 15: 
    case 16: proc(8, 12, 17); break;
    case 17: collect(4, 18, 13,14,15,16); break;
    case 18: disperse (2, 17, 19,20); break;
    case 19: 
    case 20: proc(2, 18, 21); break;
    case 21: collect(2, 22, 19,20); break;
    case 22: proc(2, 21, 23); break;
    case 23: out(no, 22); break;
  }

  MPI_Finalize();
  return 0;

}
