/************************************************
CS789 - Parallel Programming 

	Parallel implementation.

@Author: Cabel Dhoj Shrestha
Nov 2015
**************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <sys/time.h>

#define assert__(x) for ( ; !(x) ; assert(x) )

int *A, *B, *C, *T;/* Whole matrices */
int *A_hat, *B_hat, *C_hat;/* Sub-matrices */
int M, i, j, n, rank, size, N, Nsqrt;
int src, last, above, below, next, prev;
MPI_Status status;

void printTime(int rank, char *task, struct timeval t1, struct timeval t2) {

  long l;
  long secs, usecs;

  l = t2.tv_sec*1000000+t2.tv_usec-(t1.tv_sec*1000000+t1.tv_usec);

  secs = l/1000000;
  usecs = l%1000000;

  printf(">>%d:%s \t %ld.%ld\n", rank, task, secs,usecs);
}

void pretty_matrix(int *matrix, char *name, int dim) {
	printf("\n%d:Matrix %s=\n", rank, name);
	for (i=0; i<dim; i++) {
		if (i==0)
			printf(" _%*s_\n", dim*2-1," ");
		for (j=0; j<dim; j++) {
			if(j==0)
				printf("| ");
			printf("%d", matrix[j+i*dim]);
			if (j!=dim-1)
				printf(" ");
			if(j==dim-1)
				printf(" |");
		}
		printf("\n");
		if (i==dim-1)
			printf(" -%*s-\n", dim*2-1," ");
	}
}

void ReadFile(char *filename) /* Read matrices A and B */
{

	FILE *f;
	f=fopen(filename,"r"); 
	assert(f);
	for (i=0; i<M; i++)
		for (j=0; j<M; j++)
			fscanf(f,"%d ",&A[i*M+j]);

	pretty_matrix(A, "A", M);

	for (i=0; i<M; i++)
		for (j=0; j<M; j++)
			fscanf(f,"%d ",&B[i*M+j]);
	fclose(f);

	pretty_matrix(B, "B", M);
}

void WriteFileResult(char *filename) /* Write result matrix C */
{
	FILE *f;
	f=fopen(filename,"w"); 
	assert(f);
	for (i=0; i<M; i++)
	{
		for (j=0; j<M; j++)
			fprintf(f,"%d ",C[j+i*M]);
		fprintf(f,"\n");
	}
	fclose(f);
}

int compute_source(int r){
	int level = rank / Nsqrt;
	return (level * Nsqrt) + (int)((r + level) % Nsqrt);
}

int compute_last(int source){
	int level = rank / Nsqrt;
	return (Nsqrt * (rank/Nsqrt)) + (source + (Nsqrt-1)) % Nsqrt;
}

int compute_next(int r){
	int next;
	int level = rank / Nsqrt;
	if((rank+1) % Nsqrt == 0){
		next = Nsqrt * level;
	}else{
		next = rank + 1;
	}
	return next;
}

int compute_prev(int r){
	int prev;
	int level = rank / Nsqrt;

	if(rank % Nsqrt == 0){
		prev = (Nsqrt * level) + (Nsqrt-1);
	}else{
		prev = rank - 1;
	}
	return prev;
}

int compute_above(int r){
	int level = rank / Nsqrt;
	int above_level;

	if( (level-1) < 0 ){
		above_level = Nsqrt-1;
	}else{
		above_level = level - 1;
	}
	return (above_level * Nsqrt) + rank % Nsqrt;
}

int compute_below(int r){
	int below_level;
	int level = rank / Nsqrt;

	if( (level+1) > (Nsqrt-1) ){
		below_level = 0;
	}else{
		below_level = level + 1;
	}
	return (below_level * Nsqrt) + rank % Nsqrt;

}
void roll_B(int r) {

	above = compute_above(r);
	below = compute_below(r);

	printf("%d:Rolling above to %d\n", rank, above);
	MPI_Send(&B[0], M*M/N, MPI_INT, above, 0, MPI_COMM_WORLD); // Send up with wrap around

	printf("%d:Rolled down from %d\n", rank, below);
	MPI_Recv(&B[0], M*M/N, MPI_INT, below, 0, MPI_COMM_WORLD, &status); // Receive from below
}

void pipe_A(int r) {

	src = compute_source(r);
	last = compute_last(src);

	next = compute_next(r);
	prev = compute_prev(r);

	if (rank == src) {
		T = A; 
		printf("%d:Piping to %d\n", rank, last);
		MPI_Send(&T[0], M*M/N, MPI_INT, next, 0, MPI_COMM_WORLD ); // Send to right with wrap around

	} else if (rank != last) {
		printf("%d:Piped from %d\n", rank, src);
		MPI_Recv(&T[0], M*M/N, MPI_INT, prev, 0, MPI_COMM_WORLD, &status ); // Receive from left
		printf("%d: and Piping to %d\n", rank, last);
		MPI_Send(&T[0], M*M/N, MPI_INT, next, 0, MPI_COMM_WORLD ); // Send to right with wrap around
	} else {
		printf("%d:Piped from %d\n", rank, src);
		MPI_Recv(&T[0], M*M/N, MPI_INT, prev, 0, MPI_COMM_WORLD, &status); // Receive from left
	}
}

void slave_code() {

	A = (int *) calloc(M*M/N, sizeof(int));
	B = (int *) calloc(M*M/N, sizeof(int));
	C = (int *) calloc(M*M/N, sizeof(int));
	T = (int *) calloc(M*M/N, sizeof(int));

	MPI_Recv(A, M*M/N, MPI_INT, size-1, 0, MPI_COMM_WORLD, &status);
	MPI_Recv(B, M*M/N, MPI_INT, size-1, 0, MPI_COMM_WORLD, &status);

	/* Matrix multiplication */
	int r;
	for (r=0;r<Nsqrt;r++) {
		pipe_A(r);

		for (i=0; i<M/Nsqrt; i++) {
			for (j=0; j<M/Nsqrt; j++) {
				for(n=0; n<M/Nsqrt; n++)
					C[i*M/Nsqrt+j] = C[i*M/Nsqrt+j] + T[i*M/Nsqrt+n]*B[j+n*M/Nsqrt];
			}
		}

		roll_B(r);
	}
	// pretty_matrix(C, "C-result", M/Nsqrt);

	printf("%d:ALL ROUNDS DONE...SENDING TO MASTER\n", rank);

	MPI_Send(&C[0], M*M/N, MPI_INT, (size-1), 0, MPI_COMM_WORLD);
}

void master_code() {

	A = (int *) malloc(sizeof(int)*M*M);
	B = (int *) malloc(sizeof(int)*M*M);
	C = (int *) malloc(sizeof(int)*M*M);

	ReadFile("matrix.txt");

	A_hat = (int *) malloc(sizeof(int)*M*M/N);
	B_hat = (int *) malloc(sizeof(int)*M*M/N);
	C_hat = (int *) malloc(sizeof(int)*M*M/N);

	/*
	 * Sub-Matrix Division:   
	 *	if A is an M × M matrix, Aˆlk holds the elements Aij
	 *	where,
	 *		lM/sqrt(N) ≤ i < (l + 1)M/sqrt(N)
	 *		kM/sqrt(N) ≤ j < (k + 1)M/sqrt(N)
	 */
	int l, k;

	for (l=0; l<Nsqrt; l++) {
		for(k=0; k<Nsqrt; k++) {

			for(i=0; i<M/Nsqrt; i++) {
				for(j=0; j<M/Nsqrt; j++) {
					A_hat[i*M/Nsqrt+j] = A[i*M+j+M/Nsqrt*k+l*M*M/Nsqrt];
				}	
			}

			for(i=0; i<M/Nsqrt; i++) {
				for(j=0; j<M/Nsqrt; j++) {
					B_hat[i*M/Nsqrt+j] = B[i*M+j+M/Nsqrt*k+l*M*M/Nsqrt];
				}	
			}

			pretty_matrix(A_hat, "A_hat", M/Nsqrt);
			pretty_matrix(B_hat, "B_hat", M/Nsqrt);

			MPI_Send(&A_hat[0], M*M/N, MPI_INT, l*Nsqrt+k, 0, MPI_COMM_WORLD);
			MPI_Send(&B_hat[0], M*M/N, MPI_INT, l*Nsqrt+k, 0, MPI_COMM_WORLD);
		}
	}
	
	int c;
	for(l=0; l<Nsqrt; l++) {
		for(k=0; k<Nsqrt; k++) {
			MPI_Recv(&C_hat[0], M*M/N, MPI_INT, l*Nsqrt+k, 0, MPI_COMM_WORLD, &status );
			// pretty_matrix(C_hat, "C_hat", M/Nsqrt);

			// printf("slave:%d\n", l*Nsqrt+k);
			c = ((l*M + k)*M)/sqrt(N);
			for(i=0; i<M/Nsqrt; i++) {
				for(j=0; j<M/Nsqrt; j++) {
					C[c+M*i+j] = C_hat[i*M/Nsqrt+j];
				}
			}

		}
	}

	pretty_matrix(C, "C-Final", M);

	WriteFileResult("result-par.txt");
}

int main(int argc, char *argv[]) {

	struct timeval t1, t2;

  	/* --Start-time Block-- */
  	gettimeofday(&t1, NULL);
  	/*----------X-----------*/

	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	MPI_Comm_size (MPI_COMM_WORLD, &size);

	N=size-1;
	Nsqrt = sqrt(N);
	M = atoi(argv[1]);

	/* Check for perfect squareroot of slave processes.*/
	assert__(rintf(sqrt(N)) == sqrt(N)) {
		if (rank == 1)
			printf("Number of slave processes (N)=%d does not have perfect squareroot!\n", N);
		MPI_Finalize();
		exit(0);

	}

	if (rank == (size-1))
	{
		master_code();
	} else {
		slave_code();
	}

	MPI_Finalize();

	printf("*** End of Process %d ***\n", rank);

	/*------- End-time Block ---------*/
	gettimeofday(&t2, NULL);
	printTime(rank, "Total Execution", t1, t2);
	/*--------------X-----------------*/

	return 0;
}
