/*************************************************************/
/*                 POLAR TO CATRTESIAN                       */
/*************************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include "SVL_header.h"

void Pol2Cart(int New_Nx, int New_Ny, double **TH, double **RHO, double **Kx, double **Ky) { 
  int i, j;

  //printf("%s:%s:%d \n", __FILE__, __func__, __LINE__);

  for (i = 0; i < New_Nx; ++i) {
     for (j = 0; j < New_Ny; ++j)  {
         Kx[i][j] = RHO[i][j] * cos(TH[i][j]);
         Ky[i][j] = RHO[i][j] * sin(TH[i][j]);
     }    
  }
}  /*     END POL2CART     */
