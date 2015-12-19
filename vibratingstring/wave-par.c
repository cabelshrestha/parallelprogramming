/************************************************
CS789 - Parallel Programming 

  Parallel implementation of Vibrating String.

@Author: Cabel Dhoj Shrestha
Aug 2015
**************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <sys/time.h>

const double USEC = 1000000;

void printTime(int rank, char *task, int size, struct timeval t1, struct timeval t2) {

  long l;
  long secs, usecs;

  l = t2.tv_sec*1000000+t2.tv_usec-(t1.tv_sec*1000000+t1.tv_usec);

  secs = l/1000000;
  usecs = l%1000000;

  printf(">>%d:%s \t %ld.%ld\n", rank, task, secs,usecs);
}

int main(int argc, char *argv[]) {

  struct timeval tcomp1, tcomp2, tcomp;
  struct timeval tcomm1, tcomm2, tcomm = {0};
  struct timeval tw1, tw2, tw = {0};
  struct timeval tres;
  struct timeval t1, t2;

  /* --Start-time Block-- */
  gettimeofday(&t1, NULL);
  /*----------X-----------*/

  /* tcomp1 */
  gettimeofday(&tcomp1, NULL);

  if (argc != 6) {
    printf("Usage: wave <l> <nb> <n> <steps> <dt>\n");
    printf("\t l \t = Total length of the string (x-axis).\n");
    printf("\t nb \t = Number of half sine waves.\n");
    printf("\t n \t = Number of nodes (number of discrete points on the x-axis between 0 and l.)\n");
    printf("\t steps \t = Number of steps.\n");
    printf("\t dt \t = Size of each step.\n%d args supplied\n",argc);
    exit(0);
  }

  int n;
  int nb;
  long double l;
  long double dt;
  int steps;     
  FILE *y_file;
  long double pi, tau, dx;
  int i,s;
  int rank, size, received_len, block, offset; 
  int worker, master = 0;
  long double *y, *yold, *ynew, *temp;

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

  /*====================== MPI code begins from here. ======================*/
  MPI_Status status;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  /* tcomp2 */
  gettimeofday(&tcomp2, NULL);

  /* tcomp = tcomp2 - tcomp1 */
  timersub(&tcomp2, &tcomp1, &tcomp);

  /*================================ MASTER ================================*/
  if (rank == master) {
    /* tcomp1 */
    gettimeofday(&tcomp1, NULL);

    /* Allocate space for y array */
    y = (long double *) calloc(n,sizeof(long double));

    /* Hold end points of the string to x-axis */
    y[0] = 0.0;
    y[n-1] = 0.0;

    /* tcomp2 */
    gettimeofday(&tcomp2, NULL);

    /* tcomp += tcomp2 - tcomp1 */
    timersub(&tcomp2, &tcomp1, &tres);
    timeradd(&tcomp, &tres, &tcomp);

    /* Receive results from workers */
    for (worker = 1; worker < size; worker++) {
      /*
        Since, each worker is sending different length of 'y' array,
        that starts at different offset, use MPI_Probe tag and MPI_Get_count
        to find out those values dynamically instead of recalculating.
      */
      MPI_Probe(worker, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

      offset = status.MPI_TAG;
      MPI_Get_count(&status, MPI_INT, &received_len);

      /* tcomm1 */
      gettimeofday(&tcomm1, NULL);
      
      MPI_Recv(&y[offset], received_len, MPI_LONG_DOUBLE, worker, offset, MPI_COMM_WORLD, &status);
      /* tcomm2 */
      gettimeofday(&tcomm2, NULL);

      /* tcomm += tcomm2  - tcomm1 */
      timersub(&tcomm2, &tcomm1, &tres);
      timeradd(&tcomm, &tres, &tcomm);
    }

    /* tw1 */
    gettimeofday(&tw1, NULL);

    /* Write result to file */
    y_file = fopen("resultPar.txt", "w" );
    if (y_file == (FILE *)NULL) { 
      printf("Could not open output file.\n");
      MPI_Finalize();
      exit(1);
    }

    for (i=0; i<n; i++) {
      fprintf(y_file,"%Lf %15.15Lf\n",( l * i )/( n - 1 ), y[i]);
    } 

    fclose(y_file);

    /* tw2 */
    gettimeofday(&tw2, NULL);

    /* tw = tw2 - tw1 */
    timersub(&tw2, &tw1, &tw);
  } 
  /*============================== End Master ==============================*/

  /*================================ WORKER ================================*/
  else {

    /* tcomp1 */
    gettimeofday(&tcomp1, NULL);

    /* calculate optimal length of y for this worker to calculate.*/
    block = n/(size-1) + (rank <= (n % (size-1)) ? 1 : 0);

    /* Calculate the offset on y array that this worker should start from. */
    offset = (rank-1) * block + 1;

    /* allocate space for arrays */
    y = (long double *) calloc(block+2,sizeof(long double)); 
    yold = (long double *) calloc(block+2,sizeof(long double));
    ynew = (long double *) calloc(block+2,sizeof(long double));

    /* initialize y*-arrays */
    for (i = 1; i <= block; i++) {
      if (((offset+i-1) == 1) || ((offset+i-1) == n)) {
        y[i] = yold[i] = ynew[i] = 0.0;
      }
      else
        y[i] = yold[i] = sin(pi*nb*(l*(offset+i-1)/(n-1))/l);  
    }

    /* tcomp2 */
    gettimeofday(&tcomp2, NULL);

    /* tcomp += tcomp2 - tcomp1 */
    timersub(&tcomp2, &tcomp1, &tres);
    timeradd(&tcomp, &tres, &tcomp);

      /* Perform calculations */
      for (s=0; s<steps; s++) {

        /* Leftmost worker. Only sends its value to right worker. */
        if (rank == 1) {
          /* Send right point to 2 */
          MPI_Send(&y[block], 1, MPI_LONG_DOUBLE, rank+1, 0, MPI_COMM_WORLD);


          /* tcomm1 */
          gettimeofday(&tcomm1, NULL);

          /* Receive left point from 2 */
          MPI_Recv(&y[block+1], 1, MPI_LONG_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &status);

          /* tcomm2 */
          gettimeofday(&tcomm2, NULL);

          /* tcomm += tcomm2  - tcomm1 */
          timersub(&tcomm2, &tcomm1, &tres);
          timeradd(&tcomm, &tres, &tcomm);

        } 

        /* Rightmost worker. Only sends its value to left worker. */
        else if ( rank == (size -1 )) {
          /* Send left point to rank-1 */
          MPI_Send(&y[1], 1, MPI_LONG_DOUBLE, rank-1, 0, MPI_COMM_WORLD);


          /* tcomm1 */
          gettimeofday(&tcomm1, NULL);

          /* Receive right point from rank-1 */
          MPI_Recv(&y[0], 1, MPI_LONG_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &status);

          /* tcomm2 */
          gettimeofday(&tcomm2, NULL);

          /* tcomm += tcomm2  - tcomm1 */
          timersub(&tcomm2, &tcomm1, &tres);
          timeradd(&tcomm, &tres, &tcomm);
        } 

        /* All workers in between leftmost and rightmost. Sends their values to both right and left direction */
        else {

          /* Send left point to rank-1 */
          MPI_Send(&y[1], 1, MPI_LONG_DOUBLE, rank-1, 0, MPI_COMM_WORLD);

          /* Send right point to rank+1 */
          MPI_Send(&y[block], 1, MPI_LONG_DOUBLE, rank+1, 0, MPI_COMM_WORLD);

          /* tcomm1 */
          gettimeofday(&tcomm1, NULL);
          
          /* Receive left point from rank+1 */
          MPI_Recv(&y[block+1], 1, MPI_LONG_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &status);

          /* tcomm2 */
          gettimeofday(&tcomm2, NULL);
          /* tcomm += tcomm2  - tcomm1 */
          timersub(&tcomm2, &tcomm1, &tres);
          timeradd(&tcomm, &tres, &tcomm);

           /* tcomm1 */
                gettimeofday(&tcomm1, NULL);

          /* Receive right point from rank-1 */
          MPI_Recv(&y[0], 1, MPI_LONG_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &status);

          /* tcomm2 */
                gettimeofday(&tcomm2, NULL);

                /* tcomm += tcomm2  - tcomm1 */
                timersub(&tcomm2, &tcomm1, &tres);
                timeradd(&tcomm, &tres, &tcomm);

        }


        /* tcomp1 */
        gettimeofday(&tcomp1, NULL);

        /* compute */
        for (i = 1; i <= block; i++) { 
            ynew[i] = 2.0*y[i]-yold[i]+tau*tau*(y[i-1]-2.0*y[i]+y[i+1]);      
        }

        temp = yold;
        yold = y;
        y = ynew;
        ynew = temp;

        /* tcomp2 */
        gettimeofday(&tcomp2, NULL);

        /* tcomp += tcomp2 - tcomp1 */
        timersub(&tcomp2, &tcomp1, &tres);
        timeradd(&tcomp, &tres, &tcomp);

      }

      /* Send results back to master */
      MPI_Send(&y[1], block, MPI_LONG_DOUBLE, 0, offset, MPI_COMM_WORLD);

    }
    /*============================== End Worker ==============================*/

    /* tcomp1 */
    gettimeofday(&tcomp1, NULL);

    MPI_Finalize();

    /* tcomp2 */
    gettimeofday(&tcomp2, NULL);

    /* tcomp += tcomp2 - tcomp1 */
    timersub(&tcomp2, &tcomp1, &tres);
    timeradd(&tcomp, &tres, &tcomp);

    printf("rank:%d \t communication:%f \t computation:%f \t I/O:%f\n", rank, tcomm.tv_sec+tcomm.tv_usec/USEC, tcomp.tv_sec+tcomp.tv_usec/USEC, tw.tv_sec+tw.tv_usec/USEC); 

    /*------- End-time Block ---------*/
    gettimeofday(&t2, NULL);
    if (rank==0)
      printTime(rank, "Total Execution", size, t1, t2);
    /*--------------X-----------------*/

    return 0;
  }