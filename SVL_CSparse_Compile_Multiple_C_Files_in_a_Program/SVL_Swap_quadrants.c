/***************************************************/
/* SWAP QUADRANTS (1 <--> 3, 2 <--> 4) DIAGONALLY  */
/***************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <complex.h>
# include "SVL_header.h"

void SVL_SWAPQUADRANTS (int Nx, int Ny, double **A_re, double **A_im) {
  int i, j;
  int N2y, N2x; 
  double tmp13, tmp24;
   
  FILE *fp5, *fp6; // open a file
  fp5 = fopen("data5", "w"); // save Real
  fp6 = fopen("data6", "w"); // save Imagynary

  printf("\nSTEP 4: SET FFTW_SWAP\n");
  //printf("%s:%s:%d \n", __FILE__, __func__, __LINE__);
  
  N2x = Nx/2;
  N2y = Ny/2;

// Real part
  for(i = 0; i < N2x; ++i) {
     for(j = 0; j < N2y; ++j) {
	tmp13 = A_re[i][j];
	A_re[i][j] = A_re[i + N2x][j + N2y] ;
	A_re[i + N2x][j + N2y] = tmp13;

	tmp24 = A_re[i + N2x][j];
	A_re[i + N2x][j] = A_re[i][j + N2y] ;
	A_re[i][j + N2y] = tmp24;
      }
  }

// Imaginary part
  for(i = 0; i < N2x; ++i) {
     for(j = 0; j < N2y; ++j) {
	tmp13 = A_im[i][j];
	A_im[i][j] = A_im[i + N2x][j + N2y] ;
	A_im[i + N2x][j + N2y] = tmp13;

	tmp24 = A_im[i + N2x][j];
	A_im[i + N2x][j] = A_im[i][j + N2y] ;
	A_im[i][j + N2y] = tmp24;
      }
  }
  printf("\n      SVL_SWAPQUADRANTS - A_re array is pass by reference from SVL_SWAPQUADRANTS\n");
  for(i = 0; i < Nx; i++) {
    for(j = 0; j < Ny; j++) {
       fprintf(fp5,"%f   ",A_re[i][j]);
       //printf(" %f  ", A_re[i][j]);
    }
    fprintf(fp5,"\n");
  }
  
  printf("\n      SVL_SWAPQUADRANTS - A_im array is pass by reference from SVL_SWAPQUADRANTS\n");
  for(i = 0; i < Nx; i++) {
    for(j = 0; j < Ny; j++) {
       fprintf(fp6,"%f   ",A_im[i][j]);
       //printf(" %f   ", A_im[i][j]);
    }
    fprintf(fp6,"\n");
  }
    
  fclose(fp5);
  fclose(fp6);  
}/* END FUNCTION SVL_SWAPQUADRANTS */
