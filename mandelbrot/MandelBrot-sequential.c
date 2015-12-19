/************************************************
 CS789 - Parallel Programming 

  Sequential version of Mandelbrot program.
  
 @Author: Cabel Dhoj Shrestha
 Oct 2015
**************************************************/
 
#include <stdio.h>
#include <stdlib.h>

int cal_pixel(float c_real, float c_imag) {
  int count, max;
  float z_real, z_imag;
  float temp, lengthsq;

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
  return count;
}


int main(int argc, char *argv[]) {
  int disp_width, disp_height;
  float real_min, real_max, imag_min, imag_max, scale_real, scale_imag;
  FILE *f;
  int x,y,i;
  char str[256];
  float c_real, c_imag;

  int map[3][257];

  if (argc != 9) {
printf("Usage:\n Mandelbrot width height real-min real-max imag-min imag-max mapfile outfile\n");
exit(1);
  }

  /* Decode arguments */
  disp_width  = atoi(argv[1]);
  disp_height = atoi(argv[2]);

  real_min = atof(argv[3]);
  real_max = atof(argv[4]);
  imag_min = atof(argv[5]);
  imag_max = atof(argv[6]);
  
  /* Load the required colour map file */
  f = fopen(argv[7],"r");
  for(i=0;i<=256;i++) {
fgets(str,1000,f);
sscanf(str,"%d %d %d",&(map[0][i]),&(map[1][i]),&(map[2][i]));  
  }
  fclose(f);
  
  /* Compute scaling factors */
  scale_real = (real_max-real_min)/disp_width;
  scale_imag = (imag_max-imag_min)/disp_height;

  /*
   * Createing P6 format instead of P3 ppm file to speed up.
   */
  f = fopen(argv[8],"wb");
  fprintf(f,"P6\n%d %d\n255\n",disp_width,disp_height);

  for (y=0; y<disp_height; y++) {
for (x=0; x<disp_width; x++) {
  c_real = real_min + ((float) x * scale_real);
  c_imag = imag_min + ((float) y * scale_imag);
  i = cal_pixel(c_real, c_imag);

  static unsigned char color[3];
  color[0] = map[0][i];
  color[1] = map[1][i];
  color[2] = map[2][i];
  (void) fwrite(color, 1, 3, f);
}
  }

  fclose(f);
}
   