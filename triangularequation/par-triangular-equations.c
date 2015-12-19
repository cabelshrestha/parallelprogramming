/************************************************
CS789 - Parallel Programming 

  Parallel implementation of solving Upper
  Triangular System of Equations. 

@Author: Cabel Dhoj Shrestha
Sept 2015
**************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include <stdarg.h>
#include <mpi.h>

void printTime(int rank, char *task, struct timeval t1, struct timeval t2) {

  long l;
  long secs, usecs;

  l = t2.tv_sec*1000000+t2.tv_usec-(t1.tv_sec*1000000+t1.tv_usec);

  secs = l/1000000;
  usecs = l%1000000;

  printf(">>%d:%s \t %ld.%ld\n", rank, task, secs,usecs);
}

int main(int argc, char *argv[]) {

  struct timeval t1, t2;

  /* --Start-time Block-- */
  gettimeofday(&t1, NULL);
  /*----------X-----------*/
  
  int i, j, n=10;
  int rank, size;
  double sum;
  double x[10] = {0,0,0,0,0,0,0,0,0,0};
  double b[10] = {10.0, 15.0, 24.0, 5.0, 6.0, 14.0, 49.0, 11.0, 7.0, 24.0};
  double a[10][10] = {
                {2.0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                {-3.0, 12.0, 0, 0, 0, 0, 0, 0, 0, 0},
                {4.0, -5.0, 10.0, 0, 0, 0, 0, 0, 0, 0},
                {8.0, -17.0, 20.0, 5.0, 0, 0, 0, 0, 0, 0},
                {18.0, 5.0, 14.0, 9.0, -20, 0, 0, 0, 0, 0},
                {2.0, 7.0, -20.0, 8.0, 14.0, -3, 0, 0, 0, 0},
                {8.0, -6.0, -2.0, 12.0, -14.0, 9.0, -11, 0, 0, 0},
                {-9.0, 7.0, 9.0, -11.0, -2.0, 10.0, 5.0, 13, 0, 0},
                {1.0, 10.0, 6.0, 4.0, -19.0, 5.0, 22.0, 7.0, 9, 0},
                {-18.0, 5.0, -8.0, -3.0, 30.0, 12.0, 16.0, 3.0, 2.0, -8},
              };

MPI_Status status;
MPI_Init (&argc, &argv);
MPI_Comm_rank (MPI_COMM_WORLD, &rank);
MPI_Comm_size (MPI_COMM_WORLD, &size);

if (rank == 0)
{
  sum = 0;
  x[rank] = (b[rank] - sum) / a[rank][rank];
  MPI_Send(&x[rank], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
  printf("x[%d] = %f\n", rank, x[rank]);
}

else if (rank == size-1)
{
  sum = 0;
  for (j = 0; j < rank; j++) 
  {
      MPI_Recv(&x[j], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &status);
      sum = sum + a[rank][j] * x[j];
  }
  x[rank] = (b[rank] - sum) / a[rank][rank];
  printf("x[%d] = %f\n", rank, x[rank]);
}

else
{
  sum = 0;
  for (j = 0; j < rank; j++) 
  {
      MPI_Recv(&x[j], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &status);
      MPI_Send(&x[j], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
      sum = sum + a[rank][j] * x[j];
  }
  x[rank] = (b[rank] - sum) / a[rank][rank];
  MPI_Send(&x[rank], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
  printf("x[%d] = %f\n", rank, x[rank]);
}

MPI_Finalize();

/*------- End-time Block ---------*/
gettimeofday(&t2, NULL);
printTime(rank, "Total Execution", t1, t2);
/*--------------X-----------------*/

return 0;

}
