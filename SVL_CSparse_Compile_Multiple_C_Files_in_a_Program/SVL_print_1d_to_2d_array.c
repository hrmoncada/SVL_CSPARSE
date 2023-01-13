/*************************************************************/
/*      INPUT 1D ARRAY (==> PRINT ==>) OUTPUT 2D ARRAY       */
/*************************************************************/
# include <stdio.h>
# include <stdlib.h>
# include "SVL_header.h"

void PRINT_1D_ARRAY_2D( char* desc, int m, int n, double* a, int lda ) {
        int i, j; //LDA = first dimension in array, number of rows, A(M,N) then LDA=M

        //printf("%s:%s:%d \n", __FILE__, __func__, __LINE__);

        printf( "\n %s\n", desc);
        for( i = 0; i < m; i++ ) {
            for( j = 0; j < n; j++ ) {
               printf("%2.1f  ",a[i+j*lda] );
            }
            printf( "\n" );
        }
}/* END PRINT_1D_ARRAY_2D  */
