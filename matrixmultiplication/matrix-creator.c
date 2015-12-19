/************************************************
CS789 - Parallel Programming 

	Matrix generator to use with multiplication 
	program.

@Author: Cabel Dhoj Shrestha
Nov 2015
**************************************************/

#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[ ]){

	int s = atoi(argv[1]);
	FILE *f = fopen("matrix.txt","w");
	int i,j;
	srand(15);

	printf("Writing the first matrix.\n");

	for (i=0;i<s;i++) {
		for (j=0;j<s;j++)
			fprintf(f,"%d ",rand() % 3);

		// fprintf(f, "%s\n", " ");
	}

	printf("Writing the second matrix.\n");
	for (i=0;i<s;i++) {
		for (j=0;j<s;j++)
			fprintf(f,"%d ",rand() % 3);
		// fprintf(f, "%s\n", " ");
	}

	fclose(f);

	printf("Done!\n");
}
