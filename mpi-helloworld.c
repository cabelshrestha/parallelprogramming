#include <stdio.h>
#include <mpi.h>

int main(int argc, char* argv[]) {

	int rank, size;

	/* Initializing MPI. */
	MPI_Init(&argc, &argv);
	/* Getting count of total processes. */
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	/* Getting rank of the process. */
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int number;
	if (rank == 0) { //Master
		for(int i = 1; i < size; i++) {
		    number = -1;
		    MPI_Send(&number, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		    printf("back from %d process.\n",i);
		} 
	} else { //Slaves
	    MPI_Irecv(&number, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    printf("Process %d received number %d from process 0\n", rank, number);
	}

	/* Closing MPI. */
	MPI_Finalize();

	return 0;

}
