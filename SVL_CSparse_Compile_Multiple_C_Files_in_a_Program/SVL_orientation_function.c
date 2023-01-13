/***************************************************/
/*          SVL ORIENTATION FUNCTION               */
/***************************************************/
# include <stdio.h>
# include <stdlib.h>
# include "SVL_header.h"

void SVL_ORIENTATION_FUNCTION (double pi, int NN, int NM, int Lx, int Ly, int NPx, int NPy, int New_Nx, int New_Ny, double New_dx, double New_dy, double **RSQ, double **THETA) {
  
  // Calculate grating vectores 
  int i, j;
  double xa, ya, mean_xa, mean_ya;
  double X[New_Nx][New_Ny], Y[New_Nx][New_Ny];

  printf("\nSTEP 7: RSO and  THETA\n");
  //printf("%s:%s:%d \n", __FILE__, __func__, __LINE__);

  FILE  *fp9, *fp10; // open a file
  fp9 = fopen("data9", "w"); // save 
  fp10 = fopen("data10", "w"); // save
  
  for (i = 0; i < New_Nx; ++i) xa = xa + i * New_dy; 
  for (i = 0; i < New_Ny; ++i) ya = ya + i * New_dy; 

  mean_xa = (double) xa/New_Nx;
  mean_ya = (double) ya/New_Ny;
 
// meshgrid
  for (i = 0; i < New_Nx; ++i){  
      for (j = 0; j < New_Ny; ++j){
	  X[i][j] = i * (New_dx);
          THETA[i][j] = 0.0;
      }
   }

 for (j = 0; j < New_Ny; ++j){
     for (i = 0; i < New_Nx; ++i){   
          Y[i][j] = i * (New_dy);
      }
 }

 for (j = 0; j < New_Ny; ++j){
     for (i = 0; i < New_Nx; ++i){   
          RSQ[i][j] = ( X[j][i] - mean_xa ) *  ( X[j][i] - mean_xa )  +  ( Y[i][j] - mean_ya ) *  ( Y[i][j] - mean_ya );
          if ( RSQ[i][j] < 10) {
              THETA[i][j] = (double)(pi/4);//(double)(pi/4)
          }
          fprintf(fp9, " %g  " , RSQ[i][j]);
          fprintf(fp10, " %g  " , THETA[i][j]); 
      }
      fprintf(fp9, "\n");
      fprintf(fp10, "\n");
   }

fclose(fp9);
fclose(fp10);

} /*    END SVL ORIENTATION FUNCTION    */
