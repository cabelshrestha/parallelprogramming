/************************************************
CS789 - Parallel Programming 

  Sequential implementation of Vibrating String.

@Author: Cabel Dhoj Shrestha
Aug 2015
**************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void printTime(char *task, struct timeval t1, struct timeval t2) {

  long l;
  long secs, usecs;

  l = t2.tv_sec*1000000+t2.tv_usec-(t1.tv_sec*1000000+t1.tv_usec);

  secs = l/1000000;
  usecs = l%1000000;

  printf(">>%s \t %ld.%ld\n", task, secs,usecs);
}

int main(int argc, char *argv[]) {
  int n;
  int nb;
  long double l;
  long double dt;
  int steps;     
  FILE *y_file;
  long double *x, *y, *yold, *ynew;
  long double pi, tau, dx;
  int i,s;
  struct timeval t1, t2;

  /* --Start-time Block-- */
  gettimeofday(&t1, NULL);
  /*----------X-----------*/

  if (argc != 6) {
    printf("Usage: wave <l> <nb> <n> <steps> <dt>\n");
    printf("\t l \t = Total length of the string (x-axis).\n");
    printf("\t nb \t = Number of half sine waves.\n");
    printf("\t n \t = Number of nodes (number of discrete points on the x-axis between 0 and l.)\n");
    printf("\t steps \t = Number of steps.\n");
    printf("\t dt \t = Size of each step.\n%d args supplied\n",argc);
    exit(0);
  }

  /* Get variables from commando line */
  l     = atof(argv[1]);
  nb    = atoi(argv[2]);
  n     = atoi(argv[3]);
  steps = atoi(argv[4]);
  dt    = atof(argv[5]);
 
  /* set variables */
  pi = 4.0 * atan(1);
  dx = l/(n-1);
  tau = 2.0*l*dt/nb/dx;

  /* allocate space for arrays */
  x = (long double *) calloc(n,sizeof(long double));
  y = (long double *) calloc(n,sizeof(long double));
  yold = (long double *) calloc(n,sizeof(long double));
  ynew = (long double *) calloc(n,sizeof(long double));

  /* initialize x-array */
  for (i=0; i<n; i++) 
    x[i] = l*i/(n-1);

  /* initialize y*-arrays */
  for (i=0; i<n; i++) {
    if ((i==0) || (i == n-1)) 
      y[i] = yold[i] = ynew[i] = 0.0;
    else
      y[i] = yold[i] = sin(pi*nb*x[i]/l);  
  }
   
  /* Perform calculations */
  for (s=0; s<steps; s++) {
    for (i=1; i<n-1; i++) { 
      ynew[i] = 2.0*y[i]-yold[i]+tau*tau*(y[i-1]-2.0*y[i]+y[i+1]);      
    }
    for (i=1; i<n-1; i++) {
      yold[i] = y[i];
      y[i] = ynew[i];
    }
  }

  /* Write the result to a file*/
  y_file = fopen("resultSeq.txt", "w" );
  if (y_file == (FILE *)NULL) { 
    printf("Could not open output file.\n");
    exit(1);
  }
  for (i=0; i<n; i++) 
    fprintf(y_file,"%Lf %15.15Lf\n",( l * i )/( n - 1 ), y[i]);
  fclose(y_file);  

  /*------- End-time Block ---------*/
    gettimeofday(&t2, NULL);
    printTime("Total Execution", t1, t2);
    /*--------------X-----------------*/

  return 0;
}