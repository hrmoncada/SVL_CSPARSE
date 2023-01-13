/*************************************************************/
/*                   CARTESIAN TO POLAR                      */
/*************************************************************/ 
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include "SVL_header.h"

void Cart2Pol(int New_Nx, int New_Ny, double **Kx, double **Ky, double **TH ,double **RHO){ 
  int i, j;

  //printf("%s:%s:%d \n", __FILE__, __func__, __LINE__);

  for (i = 0; i < New_Nx; ++i) {
     for (j = 0; j < New_Ny; ++j)  {
         RHO[i][j] = sqrt(Kx[i][j] * Kx[i][j] + Ky[i][j] * Ky[i][j] );
         TH[i][j] =  atan2(Ky[i][j], Kx[i][j]);
     }    
  }
 }  /*   END CART2POL     */
