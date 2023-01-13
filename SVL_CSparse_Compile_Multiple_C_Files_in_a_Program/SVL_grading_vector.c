/***************************************************/
/*          SVL_GRADING_VECTOR                     */
/***************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include "SVL_header.h"

void SVL_GRADING_VECTOR (double pi, int NN, int NM, int Lx, int Ly, double *KX, double *KY)  {
  int i, j;

  printf("\nSTEP 6: (KX, KY) Wavevector cartician coord, writing KX(:) and KY(:) in a column array fashion\n");
  //printf("%s:%s:%d \n", __FILE__, __func__, __LINE__);

/* Writting Grading vectors in Column array fashion, KX(:) and KY(:) */
  for (j = 0; j < NN; ++j) {
      for (i = 0; i < NM; ++i) {
            KX[i*NM + j]  =  (2*pi/Lx)*(j - NM/2); 
	    KY[i*NM + j]  =  (2*pi/Ly)*(i - NN/2);         
        } 
    }

}/*   END GRADING VECTOR     */

