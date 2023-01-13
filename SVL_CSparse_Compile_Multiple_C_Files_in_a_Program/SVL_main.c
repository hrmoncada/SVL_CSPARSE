/***********************************************************/
/*         SPACIAL VARIANT LATTICE                         */
/***********************************************************
  Libraries:
    -lfftw3: fast fourier on the west
    -lm: math 
    -lcsparse : sparse matrix package

1. Compile:
 >> gcc -Wall -o out (write here SVL codes) -lcsparse -lfftw3 -lm 
					  
 Execute: 
 >> ./out
 >> octave plot_data.m

2. MAKE FILE: ( will run just the c file)
>> make
>> octave plot_data.m

3. RUN BASH: (compile and run everthing at ones)
   >>./SVL.sh
 or 
   >> bash SVL.sh

BASH FILE: SVL.sh
#!/bin/sh
make
./out
octave plot_data.m
**********************************************************
**********************************************************/

# include <stdio.h>
# include <stdlib.h>
//# include <fftw3.h>
# include <math.h>
# include <complex.h>
# include "cs.h"
# include <time.h>
# include "SVL_header.h"

# define MAX(a, b) (a > b ? a : b)
# define MIN(a, b) (a < b ? a : b)

/*************************************************************/
/*              START MAIN PROGRAM Program                  */
/*************************************************************/
int main() {

 clock_t tic = clock(); /* Star-t Elapse time*/

 printf("---------------------------------------\n");
 printf("         SPACIAL VARIANT LATTICE        \n");
 printf("---------------------------------------\n");
 
 int i, j;

/* Matrix size  (Nx, Ny) = (Rows, Columns) */
  int Nx = 16;
  int Ny = Nx;
      
/* Matrix length (x-axis, y-axis) */
  int Lx = 1;  
  int Ly = Lx;

/*************************************************************/
/*                 CALL CREAT MATRIX                         */
/*************************************************************/
/* U array is pass by reference from Creat_matrix*/
/* When an entire array is passed through a function, the size of the array (Nx, Ny) need to be pass as a separate parameter.*/

 double **U = CREAT_MATRIX (Nx, Ny); 

 //PRINT_2D_ARRAY_2D("2D array print as 2D, U ", Nx, Ny, U);

 FILE *fp1; // open a file
 fp1 = fopen("data1", "w"); // save mat into a file as Square Unit Cell of zeros

 for(i = 0; i < Nx; i++) {
    for(j = 0; j < Ny; j++) {
       fprintf(fp1,"%1.0f   ",U[i][j]);
    }
    fprintf(fp1,"\n");
  } 
  fclose(fp1);

/*************************************************************/
/*                 CALL SVL TRIANGLE                         */
/*************************************************************/

 SVL_TRIANGLE (Nx, Ny, U);
 
/*************************************************************/
/*                 Call SVL FFTW                             */
/*************************************************************/
 double **A_re = CREAT_MATRIX (Nx, Ny);
 double **A_im = CREAT_MATRIX (Nx, Ny);

/* A_re A_im array is pass by reference from SVL_FFTW*/
/* When an entire array is passed through a function, the size of the array (Nx, Ny) need to be pass as a separate parameter.*/

  SVL_FFTW(Nx, Ny, U, A_re, A_im); 

  free(U);
/*************************************************************/
/*                    CALL SVL IFFTW                         */
/*************************************************************/
  double **inv_A_re = CREAT_MATRIX(Nx, Ny);
  double **inv_A_im = CREAT_MATRIX(Nx, Ny);  

/* inv_A_re inv_A_im array is pass by reference from SVL_IFFTW*/
/* When an entire array is passed through a function, the size of the array (Nx, Ny) need to be pass as a separate parameter.*/

  SVL_IFFTW (Nx, Ny, A_re, A_im, inv_A_re, inv_A_im);

  free(inv_A_re);
  free(inv_A_im);
/***********************************************************/
/* CALL SVL SWAP QUADRANTS (1 <--> 3, 2 <--> 4) DIAGONALLY */
/***********************************************************/
/* A_re A_im array is pass by reference from SVL_SWAPQUADRANTS */
/* When an entire array is passed through a function, the size of the array (Nx, Ny) need to be pass as a separate parameter.*/
  
  SVL_SWAPQUADRANTS (Nx, Ny, A_re, A_im); 
    
/*************************************************************/
/*               TRUCATE THE SPATIAL HARMONICS               */
/*************************************************************/
/* Spatial Harmonics - Number of Spatial Harmonics to take into account for 2D unit cell*/
  int NM = 11; // x axis
  int NN = NM; // y axis

/* Spatial Harmonic - Truncated Array*/
  double complex *AMN  =  malloc(NM * NN * sizeof(double complex)); 

  SVL_TRUNCATE_FFTW_SPATIAL_HARMONIC (Nx, Ny, NM, NN, A_re, A_im, AMN); 

/* Free Memory */
  free(A_re);
  free(A_im); 
;
/*************************************************************/
/*                        GRADING VECTOR                     */
/*************************************************************/    
  double pi = 4.0 * atan(1.0);
  double *KX   =  malloc(NM * NN * sizeof(double));
  double *KY   =  malloc(NM * NN * sizeof(double));

  SVL_GRADING_VECTOR (pi, NN, NM, Lx, Ly, KX, KY); 
   
/*************************************************************/
/*                    DEFINE SPATIAL VARIANCE                */
/*************************************************************/
/* Lattice parameters */
  int NPx = 11;  /* Number of Unit cell to be consider*/
  int NPy = NPx;
  
/* New Matrix size */ 
  int New_Nx = 10 * NPx; //update row size matrix Nx
  int New_Ny = 10 * NPy; //update column size matrix Ny
  
  double New_dx = (double) NPx * Lx/New_Nx; //update step size dx
  double New_dy = (double) NPy * Ly/New_Ny; //update step size dy
  
  double **THETA = CREAT_MATRIX (New_Nx, New_Ny);
  double **RSQ   = CREAT_MATRIX (New_Nx, New_Ny);
  
  SVL_ORIENTATION_FUNCTION (pi, NN, NM, Lx, Ly, NPx, NPy, New_Nx, New_Ny, New_dx, New_dy, RSQ, THETA);

/*************************************************************/
/*              FDDER => DX DY DXX  DYY operator             */
/*************************************************************/
/*************************************************************/
/*           HANDLE INPUT AND OUTPUT ARGUMENTS               */
/*              M_old = 16 x 16 = Nx * Ny,                   */
/*           M_new = 110 x 110 = New_Nx * New_Ny             */
/*************************************************************/

/************************************************************************/
/*   D[m][n] is a 2D operator of dimension (2*M x M)                    */ 
/*                                                                      */
/*   m     : (int) number of rows of the operator		        */
/*   n     : (int) number of columns of the operator                    */
/*   nzmax : (int) maximum number of entries 				*/
/*   values : it is true or false   (0 or 1)     		        */
/*   triplet : it is true or false. (0 or 1) 			        */
/*                                                                      */			
/*   *p : (int) column pointers (size n+1) or col indices (size nzmax) 	*/	
/*   *i : (int) row indices, size nzmax 				*/	
/*   *j : (int) 							*/
/*   *x : (double) numerical values, size nzmax 			*/
/*   nz : (int) # of entries in triplet matrix, -1 for compressed-col	*/
/*									*/
/* allocate a sparse matrix (triplet form or compressed-column form)    */
/* add an entry to a triplet matrix; return 1 if ok, 0 otherwise        */
/*									*/
/* cs_spalloc : Creates an m-by-n sparse matrix that can hold up to     */
/*              nzmax entries. Numerical values are allocated if        */
/*              "values" is true. A triplet or compressed-column matrix */ 
/*              is allocated depending on whether "triplet" is true or  */ 
/*              false. 	         				        */
/*                                                                      */	
/*									*/
/*   cs *cs_spalloc (csi m, csi n, csi nzmax, csi values, csi triplet)  */
/*       cs_spalloc (int m, int n, int nzmax, int values, int triplet)	*/
/************************************************************************/
  int NS[2] = {New_Nx, New_Ny}; 
  int BC[2]  = {1, 1};
  int M = New_Nx * New_Ny;          // COMPUTE SIZE OF MATRICES
  double RES[2] = {New_dx, New_dy};

/* D = [DX; DY] OPERATOR IS BUILD*/
   struct cs_sparse *D   = cs_spalloc(2*M, M, 4*M, 1, 1);

   SVL_FDDER(M, NS, RES, BC, D);

/*************************************************************/
/*                    MAIN LOOP SPATIAL VARIANCE             */
/*************************************************************/
 
 SVL_LOOP (M, New_Nx, New_Ny, NM, NN, AMN, KY, KX, THETA, D);

/*************************************************************/
/*                     END MAIN LOOP                         */
/*************************************************************/ 
 
  clock_t toc = clock(); /*End Elapsed time */
  printf("\n   Elapsed: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
 
  printf("\n=============THE END===========\n");
  return 0;
}
/*************************************************************/
/*                      End main program                     */
/*************************************************************/






