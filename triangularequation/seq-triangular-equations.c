/************************************************
CS789 - Parallel Programming 

  Sequential implementation of solving Upper
  Triangular System of Equations. 

@Author: Cabel Dhoj Shrestha
Sept 2015
**************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include <stdarg.h>

void printTime(char *task, struct timeval t1, struct timeval t2) {

  long l;
  long secs, usecs;

  l = t2.tv_sec*1000000+t2.tv_usec-(t1.tv_sec*1000000+t1.tv_usec);

  secs = l/1000000;
  usecs = l%1000000;

  printf(">>%s \t %ld.%ld\n", task, secs,usecs);
}

int main(int argc, char *argv[]) {

	struct timeval t1, t2;
  
  	/* --Start-time Block-- */
  	gettimeofday(&t1, NULL);
  	/*----------X-----------*/

	int i, j, n=10;
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

	x[0] = b[0]/a[0][0]; 
	for (i = 0; i < n; i++) {
	  sum = 0;
	  for (j=0; j < i; j++) 
	    sum = sum + a[i][j] * x[j]; 
	  x[i] = (b[i] - sum) / a[i][i];
	}

	for (i = 0; i < n; i++) {
	   printf("x[%d] = %f\n", i, x[i]);
	}


	/*------- End-time Block ---------*/
	gettimeofday(&t2, NULL);
	printTime("Total Execution", t1, t2);
	/*--------------X-----------------*/

	return 0;
}
