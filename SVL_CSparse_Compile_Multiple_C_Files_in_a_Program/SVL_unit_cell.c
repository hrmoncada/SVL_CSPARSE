/*************************************************************/
/*                  ADD TRIANGLE                             */
/* INSERTING AN EQUILATERAL TRIANGLE INTO A SQUARE UNIT CELL */
/*************************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include "SVL_header.h"

void SVL_TRIANGLE (int Nx, int Ny, double **U){
 int nx1, nx2, ny1, ny2;
 int nx, ny, nxa, nxb;
 int i, j;
 double f;
 FILE *fp2; // open a file
 fp2 = fopen("data2", "w"); // save triangle into a Square Unit Cell  

 printf("\nSTEP 1: FILL U WITH ONES TO BUILD TRIANGLE-DEVICE\n");

 //printf("%s:%s:%d \n", __FILE__, __func__, __LINE__);
 
/* Compute the position indices */
  nx1 = round(0.1*Nx);
  nx2 = round(0.9*Nx);
  ny1 = round(0.1*Ny);
  ny2 = round(0.9*Ny);

  for (ny = ny1; ny < ny2; ++ny) {
      f  =  (ny - ny1)/(double)(ny2 - ny1); // f must be declare as double
      nx  = round(f*(nx2 - nx1 + 1)); // we can not round
      nxa = round((Nx - nx)/2 + 1); // me
      //nxa = round((Nx - nx)/2); // raymond
      nxb = nxa + nx - 1; 
      for (nx = nxa-1; nx <= nxb; ++nx) {
	  U[nx][ny] = 1;
      }	
  }  
  
    //printf("SVL_TRIANGLE - Add triangle inside U \n");
  for(i = 0; i < Nx; i++) {
    for(j = 0; j < Ny; j++) {
      fprintf(fp2,"  %1.0f  ", U[i][j]);  // Triangle figure 
      //printf("%1.0f   ",U[i][j]);
    }  
    fprintf(fp2, "\n");	
    //printf("\n");
  } 
  
  fclose(fp2);
  
} /* END FUNCTION SVL_TRIANGLE  */
