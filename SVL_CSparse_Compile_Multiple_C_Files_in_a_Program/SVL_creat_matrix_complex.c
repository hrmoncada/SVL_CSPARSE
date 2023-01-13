/*************************************************************/
/*                 Call Creat Matrix Complex                 */
/*************************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <complex.h>
# include "SVL_header.h"

double **CREAT_MATRIX_COMPLEX(int rows, int cols) {
 int j; 
 double complex **mat; /* Defining a temporal "mat" array , otherwise would have to use (*memory) everywhere U is used (yuck) */

//printf("%s:%s:%d \n", __FILE__, __func__, __LINE__);

 /* Each row should only contain double*, not double**, because each row will be an array of double */
 mat = (double complex **) malloc(rows * sizeof(double complex*)); // create Nx-row temporal pointers array 

 if (mat == NULL) {
    printf("Failure to allocate room for row pointers.\n ");
    exit(0);
 }

/* Had an error here.  Alloced rows above so iterate through rows not cols here */
   for (j = 0; j < cols; j++) {
/* Allocate array, store pointer  */
     mat[j] = (double complex *) malloc(cols * sizeof(double complex)); // Create on each temporal rows we allocate Ny-columns temporal pointer
    
      if (mat[j] == NULL) {
            printf("Failure to allocate for row[%d]\n", j);
            exit(0); 
/* still a problem here, if exiting with error, should free any column mallocs that were successful. */
       } 
 } 
 
  return (mat); // return mat as U
  free(mat); // free mat No U
} /* END FUNCTION CREAT MATRIX */
