/************************************************
CS789 - Parallel Programming 

	Sequential implementation.

@Author: Cabel Dhoj Shrestha
Nov 2015
**************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <sys/time.h>


int M, i, j, n;
int *A, *B, *C;

void printTime(char *task, struct timeval t1, struct timeval t2) {

  long l;
  long secs, usecs;

  l = t2.tv_sec*1000000+t2.tv_usec-(t1.tv_sec*1000000+t1.tv_usec);

  secs = l/1000000;
  usecs = l%1000000;

  printf(">>%s \t %ld.%ld\n", task, secs,usecs);
}

void pretty_matrix(int *matrix, char *name, int dim) {
	printf("\nMatrix %s=\n", name);
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
			fscanf(f,"%d ",&A[j+i*M]);

	pretty_matrix(A, "A", M);

	for (i=0; i<M; i++)
		for (j=0; j<M; j++)
			fscanf(f,"%d ",&B[j+i*M]);
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

int main(int argc, char *argv[ ]){

	struct timeval t1, t2;

  	/* --Start-time Block-- */
  	gettimeofday(&t1, NULL);
  	/*----------X-----------*/

	M = atoi(argv[1]);
	A = (int*) malloc(sizeof(int)*M*M);
	B = (int*) malloc(sizeof(int)*M*M);
	C = (int*) malloc(sizeof(int)*M*M);	

	ReadFile("matrix.txt");
	for (i=0; i<M; i++) {
		for (j=0; j<M; j++) {
			C[j+i*M] = 0;
			for(n=0; n<M; n++)
				C[j+i*M] = C[j+i*M] + A[i*M+n]*B[j+n*M];
		}
	}

	pretty_matrix(C, "C-Final", M);

	WriteFileResult("result-seq.txt");

	/*------- End-time Block ---------*/
	gettimeofday(&t2, NULL);
	printTime("Total Execution", t1, t2);
	/*--------------X-----------------*/

	return 0;

}
