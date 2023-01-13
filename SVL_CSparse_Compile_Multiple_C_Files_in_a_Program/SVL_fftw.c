/***************************************************/
/*                  FFTW                           */
/*     FAST FOURIER TRANSFORM ON THE WEST          */
/***************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <fftw3.h>
# include <math.h>
# include <complex.h>
# include "SVL_header.h"

void SVL_FFTW(int Nx, int Ny, double **U, double **A_re, double **A_im) {
  int i, j;   
  FILE *fp3, *fp4; // open a file
  fp3 = fopen("data3", "w"); // save Real
  fp4 = fopen("data4", "w"); // save Imagynary

  printf("\nSTEP 2: SET FFTW DATA\n\n");
  //printf("%s:%s:%d \n", __FILE__, __func__, __LINE__);

  fftw_plan My_plan_fft;
// pointers to input and output data 
  fftw_complex  *U_in , *U_out; 
  
// Allocate input arrays
   U_in  = fftw_malloc(sizeof(fftw_complex) * Nx * Ny);
   
// Allocate output arrays
   U_out = fftw_malloc(sizeof(fftw_complex) * Nx * Ny);

// input data fft
 for (i = 0; i < Nx; ++i) {
     for (j = 0; j < Ny; ++j) {
      U_in[i*Ny + j][0] = U[i][j]; // real part
      U_in[i*Ny + j][1] = 0.0;     // imaginary part
      } 
 }

// Create plans for fftw 2D. Here, we set the type of transformation we want, i.e. 1D, 2D or 3D and fft(FORWARD) or ifft(BACKWARD)
   My_plan_fft = fftw_plan_dft_2d(Nx, Ny, U_in, U_out, FFTW_FORWARD, FFTW_ESTIMATE);

// Execute fftw plan, output file U_out 
   fftw_execute(My_plan_fft);

// Normalize values 
   double normalization = (double) (Nx * Ny); 
   
   for (i = 0; i < Nx * Ny; i++) {
       U_out[i][0] = U_out[i][0] / normalization;
       U_out[i][1] = U_out[i][1] / normalization; 
   }

// output data fft  
   for (i = 0; i < Nx; ++i) {
       for (j = 0; j < Ny; ++j) {
           A_re[i][j] =  U_out[i*Ny + j][0];
           A_im[i][j] =  U_out[i*Ny + j][1];
       }
   }

// Destroy plan 
  fftw_destroy_plan(My_plan_fft);  
  fftw_free(U_in);		   //Free memory
  fftw_free(U_out);		   //Free memory
  
  printf("      SVL_FFTW - A_re array is pass by reference from SVL_FFTW\n"); 
  for(i = 0; i < Nx; i++) {
    for(j = 0; j < Ny; j++) {
       fprintf(fp3,"%f   ",A_re[i][j]);
       //printf("%f   ",A_re[i][j]);
    }
    fprintf(fp3,"\n");
    //printf("\n");
  }
  
  printf("\n      SVL_FFTW - A_im array is pass by reference from SVL_FFTW\n");
  for(i = 0; i < Nx; i++) {
    for(j = 0; j < Ny; j++) {
       fprintf(fp4,"%f   ",A_im[i][j]);
       //printf("%f   ",A_im[i][j]);
    }    
    fprintf(fp4,"\n");
    //printf("\n");
  } 

  fclose(fp3);
  fclose(fp4);

} /* END FUNCTION SL_FFTW */

