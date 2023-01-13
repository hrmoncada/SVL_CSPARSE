/***************************************************/
/*     TRUCATE THE SPATIAL HARMONICS               */
/***************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <complex.h>
# include "SVL_header.h"

void SVL_TRUNCATE_FFTW_SPATIAL_HARMONIC (int Nx, int Ny, int NM, int NN, double **A_re, double **A_im, double complex *AMN) { 
// notice the function is void as it returns nothing
// you passed in the address of the variables
// pointers are a pointer to an address
// these statements are basically saying that the address contains the following value
  int i, j;

  printf("\nSTEP 5: SET TRUNCATION ARRAY OF FFTW SPATIAL HARMONICS\n");
  //printf("%s:%s:%d \n", __FILE__, __func__, __LINE__);
  
  int m0, n0, n1, n2, m1, m2; 
                          // floor: return the last integral value less than o equal to x
  m0 = 1 + floor(Nx/2);  // x-axis middle point, ex: floor of 12.4 is 12.0
  n0 = 1 + floor(Ny/2);  // y-axis middle point, ex: floor of 12.4 is 12.0

  m1 = m0 - floor(NM/2) - 1; // ex: floor of 12.6 is 12.0
  m2 = m0 + floor(NM/2) - 1; // ex: floor of 13.1 is 13.0
  n1 = n0 - floor(NN/2) - 1; // ex: floor of -2.3 is -3.0
  n2 = n0 + floor(NN/2) - 1; // ex: floor of -3.8 is -4.0
   
  //int array[ ] = { m1, m2, n1, n2 };
  //int array_size  = 4;

  printf("\n         Truncated Harmonic Values         \n");
  printf("(m0, m1, m2, n0, n1, n2) = (%d, %d, %d, %d, %d, %d)\n", m0, m1, m2, n0, n1, n2);
  printf("\n         Vector Column AMN(:) \n");

  double complex **TA = CREAT_MATRIX_COMPLEX (NM, NN);
   
  FILE *fp7, *fp8; //open a file
  fp7 = fopen("data7", "w"); 
  fp8 = fopen("data8", "w");  
  //int n1, n2, m1, m2;  

 /* m1 = array[0]; // ex: floor of 12.6 is 12.0
  m2 = array[1]; // ex: floor of 13.1 is 13.0
  n1 = array[2]; // ex: floor of -2.3 is -3.0
  n2 = array[3]; // ex: floor of -3.8 is -4.0*/
  
// MAGNITUDE
  for (i = m1 ; i <= m2; ++i) {
      for (j = n1 ; j <= n2; ++j) {
	   //TA[i-m1][j-n1] = sqrt((A_re[i][j] * A_re[i][j]) + (A_im[i][j] * A_im[i][j]));
           TA[i-m1][j-n1] = A_re[i][j] + I*A_im[i][j];
           fprintf(fp7,"%f   ",  creal(TA[i-m1][j-n1])); 
           fprintf(fp8,"%f   ",  cimag(TA[i-m1][j-n1]));           
       }	
       fprintf(fp7,"\n");
       fprintf(fp8,"\n");
   } 
  
/* Spatial Harmonics - Rewritting TA as a Column Array, AMN = TA(:) */
  for (j = 0; j < NN; ++j) {
       for (i = 0; i < NM; ++i) {
	   AMN[i*NM + j]  =   TA[i][j]; 
      }	
   } 

   fclose(fp7); // close a file
   fclose(fp8); // close a file
   free(TA);
} //END FUNCTION TRUNCATE FFTW SPATIAL HARMONIC 
