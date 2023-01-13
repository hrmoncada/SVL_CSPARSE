/*************************************************************/
/*     INPUT 2D ARRAY (==> PRINT ==>) OUTPUT 2D ARRAY        */
/*************************************************************/
# include <stdio.h>
# include <stdlib.h>
# include "SVL_header.h"

void PRINT_2D_ARRAY_2D(char* desc, int m, int n, double **a) {    
   int i, j;

   //printf("%s:%s:%d \n", __FILE__, __func__, __LINE__);

   printf( "\n %s\n", desc );

   for (i = 0; i < m; ++i) {
       for (j = 0; j < n; ++j) { 
           printf("%2.1f  ",a[i][j]);
       }
      printf("\n");	
   }
} /* END PRINT_2D_ARRAY_2D  */ 
