/************************************************
CS789 - Parallel Programming 

	Parallel program to compute partial sum.

@Author: Cabel Dhoj Shrestha
Oct 2015
**************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <sys/time.h>
#include <stdarg.h>
#include <mpi.h>

int main(int argc, char *argv[ ]) {
	int rank, size;
	int *A;
	int N = atoi(argv[1]);

	/* MPI code begins here */
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);


	if (rank == 0) { /* Master */
		int i;
		A = (int *) calloc(N,sizeof(int));
		/* Initializing array A with 1...N values */
		for (i=1;i<=N;i++)
			A[i-1] = i;

		/* Send A to root node (Rank 1) of the distribution binary 
		 * tree. It will take the first value of the array and pass
		 * the remaining halfs to left (first half) and right (second half).
		 */
		MPI_Send(&A[0], N, MPI_INT, 1, 0, MPI_COMM_WORLD);

	} else { /* Slave */

		/* lengths */
		int rl, hl;
		/* nodes */
		int ln, rn;
		/* node values */
		int nv, nps, fps;
		/* others */
		int i, offset, temp, steps;

		MPI_Status status;
		MPI_Probe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
      	MPI_Get_count(&status, MPI_INT, &rl);
      	A = (int *) calloc(rl,sizeof(int));
      	MPI_Recv(&A[0], rl, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD, &status );

      	/* Take the first element as node value */
      	nv = A[0];
      	printf("%d:Assigned Value: %d\n", rank, nv);

      	/* Leaf nodes will get just 1 value each. No division needed. */
      	if (rl != 1) {
	      	/* Send halves to left and right. */
	      	hl = (rl - 1)/2;
	      	ln = rank + 1;
	      	rn = rank + hl + 1;
      	}

      	steps = (int) log2((double) (N+1));
      	nps = nv;
      	offset = 1;

      	/* Begin the partial sum calculations. */
      	for (i=0; i<steps; i++) {

      		/* Each process sends its sum to process on the right w/ offset 1, 2, 4, 8 and so until last process.*/
        	if ((rank + offset) <= N)
          		MPI_Ssend(&nps, 1, MPI_INT, rank+offset, 0, MPI_COMM_WORLD);

          	/* Each process receives from processes to its left. Any process below offset is done with its sum.*/
        	if (rank > offset) {
          		MPI_Recv(&temp, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
          		nps += temp;
        	}

        	/* Incrementing the offset for the next round.*/
        	offset *= 2;
      	}


	    if (rank == N) { /* Final node has the total partial sum. Calculation is complete, so pass the sum to all others.*/
	    	fps = nps;
	        MPI_Ssend(&fps, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
	    } else { /* Receive the final partial sum and continue distributing it. */

	    	MPI_Recv(&fps, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	        if (rl != 1) {
	        	MPI_Ssend(&fps, 1, MPI_INT, ln, 0, MPI_COMM_WORLD);
          	if (rn != N)
            	MPI_Ssend(&fps, 1, MPI_INT, rn, 0, MPI_COMM_WORLD);
	        }
	    }

    	printf ("%d:node_partial_sum=%d, final_partial_sum=%d\n", rank, nps, fps);
	}

	MPI_Finalize();
	return 0;
}