/***************************************************/
/*                  IFFTW                          */
/***************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <fftw3.h>
# include <math.h>
# include <complex.h>
# include "SVL_header.h"

void SVL_IFFTW (int Nx, int Ny, double **A_re, double **A_im, double **inv_A_re, double **inv_A_im) {
  int i, j;

  printf("%s:%s:%d \n", __FILE__, __func__, __LINE__);

  printf("\nSTEP 3: SET IFFTW DATA\n"); 
    
  fftw_plan My_plan_ifft;
  fftw_complex *inv_U_in, *inv_U_out;  // pointers to input and output data
 
// Allocate input arrays
  inv_U_in  = fftw_malloc(sizeof(fftw_complex) * Nx * Ny);
  
// Allocate output arrays
  inv_U_out = fftw_malloc(sizeof(fftw_complex) * Nx * Ny);

// Creat ifftw plan for fftw 2D. Here, we set what kind of transformation we want, i.e. 1D, 2D or 3D and fft(FORWARD) or ifft(BACKWARD)
  My_plan_ifft = fftw_plan_dft_2d(Nx, Ny, inv_U_in, inv_U_out, FFTW_BACKWARD, FFTW_ESTIMATE);

  for (i = 0; i < Nx; ++i) {
      for (j = 0; j < Ny; ++j) {
          inv_U_in[i*Ny + j][0] = A_re[i][j]; // real part
          inv_U_in[i*Ny + j][1] = A_im[i][j]; // imaginary part
      } 
  }

// Execute ifftw plan, output file inv_U_out 
   fftw_execute(My_plan_ifft);

//output data ifft
/*   for (i = 0; i < Nx; ++i) {
       for (j = 0; j < Ny; ++j) {
           inv_A_re[i][j] =  inv_U_out[i*Ny + j][0];
           inv_A_im[i][j] =  inv_U_out[i*Ny + j][1];
      //printf("{ %g, %g }  ", inv_A_re[i][j], inv_A_im[i][j]);
       }
    //printf("\n");
   }
*/
// Destroy  inv plan 
  fftw_destroy_plan(My_plan_ifft);  
  fftw_free(inv_U_in);		   //Free memory
  fftw_free(inv_U_out);		   //Free memory
} /* END FUNCTION SL_IFFTW */
