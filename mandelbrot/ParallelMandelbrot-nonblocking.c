/************************************************
 CS789 - Parallel Programming 

  A simple MPI program to create Mandelbrot picture using non blocking method. To achieve
  the non blocking nature, MPI_Iprobe is used which checks if a worker(any) is done with
  its calculation.

  The program contains one master (process 0) and multiple workers. Each worker calculates
  the amount of work the needs to be done by dividing the total height by the number of 
  workers in the program. To make it optimal, the remainder of the division is distributed
  between the workers by 1 or 0 row each.

  After each worker is done with the calculation, it sends its block of work to the master. 
  Master writes the Mandelbrot file after all workers are done.
  
@Author Cabel Dhoj Shrestha
Oct 2015
*************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <sys/time.h>

/*
 Function to calculate real time as UNIX time command does not work for parallel programs.
*/
void printTime(int rank, char *task, int size, struct timeval t1, struct timeval t2) {

  long l;
  long secs, usecs;

  l = t2.tv_sec*1000000+t2.tv_usec-(t1.tv_sec*1000000+t1.tv_usec);

  secs = l/1000000;
  usecs = l%1000000;

  printf("%d:%s, time:,%d,%ld.%ld\n", rank, task, size, secs,usecs);
}

/*
 The main method.
*/
int main(int argc, char *argv[]) {

  if (argc != 9) {
    printf("Usage:\n Mandelbrot width height real-min real-max imag-min imag-max mapfile outfile\n");
    exit(1);
  }

  struct timeval t1, t2, comm1, comm2, compu1, compu2, write1, write2;

  /*------- Start-time Block ---------*/
  gettimeofday(&t1, NULL);
  /*--------------X-----------------*/

  /*------ VARIABLES DECLARATION ------*/
  int disp_width, disp_height;
  float real_min, real_max, imag_min, imag_max, scale_real, scale_imag;
  float c_real, c_imag;
  int count, max;
  float z_real, z_imag;
  float temp, lengthsq;

  FILE *f;
  int x,y,i,worker;
  char str[256];
  int rank, size, block_height;
  int *block; 
  int received_from;
  MPI_Status recv_status;
  /*----- End of variable declarations. -----*/

  /* Decode arguments */
  disp_width  = atoi(argv[1]);
  disp_height = atoi(argv[2]);
  real_min = atof(argv[3]);
  real_max = atof(argv[4]);
  imag_min = atof(argv[5]);
  imag_max = atof(argv[6]);

  /*============== MPI code begins from here. ==============*/
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  /*================================ MASTER ================================*/
  if (rank == 0) {
    int worker_count = 0;
    int incoming=0;
    int received_len;
    int map[3][257];

    int *pic = (int *)malloc(disp_width*disp_height*sizeof(int*));
    
    /*===================Non-blocking Code ========================*/
    while(worker_count != (size-1)) {
      MPI_Status probe_status;

      /*
        Probe to seeif there is any worker done with their block and 
        is trying to send to master

        Probe time is very minimal. But since it is non-blocking, 
        there are too many outputs for this. So, ignoring time
        calculation for this.
      */
      MPI_Iprobe(MPI_ANY_SOURCE,0,MPI_COMM_WORLD,&incoming,&probe_status);

      if(incoming) {
        MPI_Get_count(&probe_status, MPI_INT, &received_len);

        block_height = received_len / disp_width;
        block = (int*)malloc(block_height * disp_width * sizeof(int*));

        /*------- Start-time Block ---------*/
        gettimeofday(&comm1, NULL);
        /*--------------X-----------------*/

        MPI_Recv(&block[0], block_height * disp_width, MPI_INT, probe_status.MPI_SOURCE, 0, MPI_COMM_WORLD, &recv_status);

        /*------- End-time Block ---------*/
        gettimeofday(&comm2, NULL);
        printTime(rank, "Communication-recv", size, comm1, comm2);
        /*--------------X-----------------*/

       /* Write the block calculation to pic array to write the Mandelbrot later.*/
        for (y=0; y<block_height; y++) {
          for (x=0; x<disp_width; x++) {
            pic[((y + ((probe_status.MPI_SOURCE-1)*block_height)) * disp_width) + x] = block[(y * disp_width) + x];
          }
        }
        worker_count++;
        incoming=0;
      }
    }

    /*------- Start-time Block ---------*/
    gettimeofday(&write1, NULL);
    /*--------------X-----------------*/

    /* Load the required colour map file */
    f = fopen(argv[7],"r");
    for(i=0;i<=256;i++) {
      fgets(str,1000,f);
      sscanf(str,"%d %d %d",&(map[0][i]),&(map[1][i]),&(map[2][i]));  
    }
    fclose(f);

    /* Creating P6 format instead of P3 ppm file to speed up. */
    f = fopen(argv[8],"wb");
    fprintf(f,"P6\n%d %d\n255\n",disp_width,disp_height);

    for (y=0; y<disp_height; y++) {
      for (x=0; x<disp_width; x++) {
        static unsigned char color[3];
        color[0] = map[0][pic[(y * disp_width) + x]];
        color[1] = map[1][pic[(y * disp_width) + x]];
        color[2] = map[2][pic[(y * disp_width) + x]];
        (void) fwrite(color, 1, 3, f);
      }
    }
    fclose(f);

    free(pic);
    /*------- End-time Block ---------*/
    gettimeofday(&write2, NULL);
    printTime(rank, "I/O", size, write1, write2);
    /*--------------X-----------------*/
  } 
  /*============================== End Master ==============================*/

  /*================================ SLAVES ================================*/
  else {
    /*------- Start-time Block ---------*/
    gettimeofday(&compu1, NULL);
    /*--------------X-----------------*/

    /* Compute scaling factors */
    scale_real = (real_max-real_min)/disp_width;
    scale_imag = (imag_max-imag_min)/disp_height;

    /* Calculate work to be done. */
    block_height = disp_height/(size - 1) + (rank <= disp_height%(size-1) ? 1 : 0);

    /* Allocate memory. */
    block = (int*)malloc(disp_width * block_height * sizeof(int*));

    /* Perform Mandelbrot calculation */
    for (y= 0; y< block_height; y++) {
      for (x=0; x<disp_width; x++) {
        c_real = real_min + (x * scale_real);
        c_imag = imag_min + ((y + ((rank-1) * block_height)) * scale_imag);

        max = 256;
        z_real = 0;
        z_imag = 0;
        count =0;
        do {
          temp = z_real * z_real - z_imag * z_imag + c_real;
          z_imag = 2 * z_real * z_imag + c_imag;
          z_real = temp;
          lengthsq = z_real * z_real + z_imag * z_imag;
          count++;
        } while ((lengthsq < 4.0) && (count < max));


        block[(y * disp_width) + x] = count;
      }
    }

    /*------- End-time Block ---------*/
    gettimeofday(&compu2, NULL);
    printTime(rank, "Computation", size, compu1, compu2);
    /*--------------X-----------------*/

    MPI_Send(&block[0], block_height * disp_width, MPI_INT, 0, 0, MPI_COMM_WORLD );
  }
  /*============================== End Slaves ==============================*/


  /*------- End-time Block ---------*/
  gettimeofday(&t2, NULL);
  printTime(rank, "Total Execution", size, t1, t2);
  /*--------------X-----------------*/

  free(block);

  MPI_Finalize();

  return 0;
}

