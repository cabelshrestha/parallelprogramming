/************************************************
 CS789 - Parallel Programming 

   A simple MPI program to create Mandelbrot picture using worker farm method.

   The program contains one master (process 0) and multiple workers. The master
   process distributes blocks of work each to the workers to begin with. 
   i.e. a chunk of height of the picture define by program argument. After each
   process is done with the calculation for its assigned block, the master sends
   another block if there is more work.

   After calculations of the whole picture is done, master writes Mandelbrot to
   a file.

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

	printf(">>%d:%s, time:,%d,%ld.%ld\n", rank, task, size, secs,usecs);
}

/*
 The main method.
*/
int main(int argc, char *argv[]) {

	if (argc != 10) {
		printf("Usage:\n Mandelbrot width height real-min real-max imag-min imag-max mapfile outfile\n");
		exit(1);
	}  

	struct timeval t1, t2, tcomm1, tcomm2, tcompu1, tcompu2, write1, write2;

	/* --Start-time Block-- */
	gettimeofday(&t1, NULL);
	/*----------X-----------*/

	/*------ VARIABLES DECLARATION ------*/
	// Variables for Mandelbrot picture calculations
	int disp_width, disp_height;
	float real_min, real_max, imag_min, imag_max, scale_real, scale_imag;
	float c_real, c_imag;
	int count, max;
	float z_real, z_imag;
	float temp, lengthsq;

	// Variables for data distribution
	int *block; 
	int block_height, offset, limit;

	// Index variables for loops
	int x,y,i,worker;

	// Variables for MPI
	MPI_Status status;
	int reply[2];
	int rank, size;

	// Variables for flags
	int more_work = 1;
	/*----- End of variable declarations. -----*/

	/* Decode arguments */
	disp_width  = atoi(argv[1]);
	disp_height = atoi(argv[2]);
	block_height = atoi(argv[3]);
	real_min = atof(argv[4]);
	real_max = atof(argv[5]);
	imag_min = atof(argv[6]);
	imag_max = atof(argv[7]);

	/*============== MPI code begins from here. ==============*/
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	/*========== MASTER =========*/
	if (rank == 0) {
		int total_sent = 0;
		char str[256];
		int map[3][257];
		FILE *f;
		int received_len, received_block, received_from;
		int active_processes = 0;
		int incoming = 0;

		int *pic = (int *)malloc(disp_width*disp_height*sizeof(int*));;

		/* Initial work distribution. 'Work' here is part of the picture height. */
		for(worker = 1; worker < size; worker++) {
			if(block_height <= disp_height) {

				if ((disp_height-total_sent) < block_height) {
					reply[0] = total_sent;
					reply[1] = disp_height;
					total_sent += disp_height - total_sent;
					more_work = 0;
				} else {
					reply[0] = (worker -1) * block_height;
					reply[1] = reply[0] + block_height;
					total_sent += block_height;
				}
			} else {
				reply[0] = 0;
				reply[1] = disp_height;
				total_sent = disp_height;
				more_work = 0;
			}

			MPI_Send(&reply[0], 2, MPI_INT, worker, 0, MPI_COMM_WORLD );
			active_processes++;

			if(!more_work) {
				more_work = 1;
				break;
			}
		}

		/* 
		   Master continues to send/recv as long as
		   there is more work to be done or there are
		   workers that haven't completed their calculation.
		*/
		while(more_work || active_processes!= 0) {

			 /*
			 	Probe to seeif there is any worker done with their block and 
			 	is trying to send to master

        		Probe time is very minimal. But since it is non-blocking, 
        		there are too many outputs for this. So, ignoring time
        		calculation for this.
      		*/
			MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&incoming,&status);

			/* Receive from worker */
			if (incoming) {
				received_from = status.MPI_SOURCE;
				offset = status.MPI_TAG;
				MPI_Get_count(&status, MPI_INT, &received_len);

				received_block = received_len / disp_width; 
				block = (int*)malloc(disp_width * received_block * sizeof(int*));

				/* --Start-time Block-- */
				gettimeofday(&tcomm1, NULL);
				/*----------X-----------*/

				MPI_Recv(&block[0], received_block * disp_width, MPI_INT, received_from, offset, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

				/*------- End-time Block ---------*/
				gettimeofday(&tcomm2, NULL);
				printTime(rank, "Communication:recv", size, tcomm1, tcomm2);
				/*--------------X-----------------*/

				active_processes--;

				/* Write the block calculation to pic array to write the Mandelbrot later.*/
				for (y=0; y<received_block; y++) {
					for (x=0; x<disp_width; x++) {
						pic[((y + offset) * disp_width) + x] = block[(y * disp_width) + x];
					}
				}

				/* Check to see if there is more work and send to worker. */
				if(total_sent == disp_height) {
					more_work = 0;
				} else if ((disp_height-total_sent) < block_height){
					reply[0] = total_sent;
					reply[1] = disp_height;
					total_sent += disp_height - total_sent;
				} else  {
					reply[0] = total_sent;
					reply[1] = total_sent + block_height;
					total_sent += block_height;
				}

				if (more_work==1) {
					/* --Start-time Block-- */
					gettimeofday(&tcomm1, NULL);
					/*----------X-----------*/

					MPI_Send(&reply[0], 2, MPI_INT, received_from, 0, MPI_COMM_WORLD );

					/*------- End-time Block ---------*/
					gettimeofday(&tcomm2, NULL);
					printTime(rank, "Communication:send-work", size, tcomm1, tcomm2);
					/*--------------X-----------------*/

					active_processes++;
				} else {
					/* 
					   If there is no more work and master has received from all workers,
					   send message to everyone telling them to terminate.
					*/
					if(active_processes==0) {
						for (worker=1; worker<size; worker++) {
							reply[0] = -1;
							reply[1] = 0;

							/* --Start-time Block-- */
							gettimeofday(&tcomm1, NULL);
							/*----------X-----------*/

							MPI_Send(&reply[0], 2, MPI_INT, worker, 0, MPI_COMM_WORLD );

							/*------- End-time Block ---------*/
							gettimeofday(&tcomm2, NULL);
							printTime(rank, "Communication:send-terminate", size, tcomm1, tcomm2);
							/*--------------X-----------------*/
						}
					}
				}
				incoming = 0;

			}
		}

		/*------- Start-time Block ---------*/
		gettimeofday(&write1, NULL);
		/*--------------X-----------------*/

		/* Load the required colour map file */
		f = fopen(argv[8],"r");
		for(i=0;i<=256;i++) {
			fgets(str,1000,f);
			sscanf(str,"%d %d %d",&(map[0][i]),&(map[1][i]),&(map[2][i]));  
		}
		fclose(f);

		/* Writing Mandelbrot to a ppm file using P6 format. */
		f = fopen(argv[9],"wb");
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
		gettimeofday(&tcompu1, NULL);
		/*--------------X-----------------*/

		/* Compute scaling factors */
		scale_real = (real_max-real_min)/disp_width;
		scale_imag = (imag_max-imag_min)/disp_height;

		while (more_work) {

			/* Receive work from master */
			MPI_Recv(&reply[0], 2, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

			offset = reply[0];
			limit = reply[1];

			/* If master sends an offset of -1, it means there is no more work to be done.*/
			if (offset == -1) {
				more_work = 0;
			} else {

				/* Perform Mandelbrot calculation */
				block = (int*)malloc(disp_width * (limit-offset) * sizeof(int*));
				for (i=0, y= offset; y< limit; i++,y++) {
					for (x=0; x<disp_width; x++) {
						c_real = real_min + (x * scale_real);
						c_imag = imag_min + (y * scale_imag);

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

						block[(i * disp_width) + x] = count;
					}
				}
				/*------- End-time Block ---------*/
				gettimeofday(&tcompu2, NULL);
				printTime(rank, "Computation", size, tcompu1, tcompu2);
				/*--------------X-----------------*/

				/*Send the calculated block to master. */
				MPI_Send(&(block[0]), (limit-offset) * disp_width, MPI_INT, 0, offset, MPI_COMM_WORLD );

			}
		}
	}
	/*============================== End Slaves ==============================*/

	/*------- End-time Block ---------*/
	gettimeofday(&t2, NULL);
	printTime(rank, "Total Execution", size, t1, t2);
	/*--------------X-----------------*/

	// free(block);

	MPI_Finalize();

	return 0;
}

