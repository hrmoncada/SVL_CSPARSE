/***********************************************************/
/*         SPACIAL VARIANT LATTICE                         */
/***********************************************************
 fftw_test_2.c
  
  Libraries:
    -lfftw3: fast fourier on the west
    -lm: math 
    -lcsparse : sparse matrix package

 Compile:
 >> gcc -Wall -o out C_Code_SVL_Mar_11_2015_CSparse.c -lcsparse -lfftw3 -lm 

						  
 Execute:
 
 >> ./out

RUN BASH: 
   >>./SVL.sh
 or 
   >> bash SVL.sh

FILE: SVL.sh
#!/bin/sh
make
./out
octave plot_data.m

MAKE FILE: ( will run just the c file)
>> make

***********************************************************
Matrix - arrangement of mathematical elements: A rectangular array of mathematical elements, e.g. the coefficients of linear equations, whose rows and columns can be combined with those of other arrays to solve problems.

Array - data structure: An arrangement of items of computerized data in tabular form for easy reference. A computer program references an item by naming the array and the position of the item in it.

MATLAB® :
MATLAB has two different types of arithmetic operations: array operations and matrix operations. You can use these arithmetic operations to perform numeric computations, for example, adding two numbers, raising the elements of an array to a given power, or multiplying two matrices.

Matrix operations follow the rules of linear algebra. By contrast, array operations execute element by element operations and support multidimensional arrays. The period character (.) distinguishes the array operations from the matrix operations. However, since the matrix and array operations are the same for addition and subtraction, the character pairs .+ and .- are unnecessary.
**********************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <fftw3.h>
# include <math.h>
# include <complex.h>
# include <limits.h>
# include "cs.h"
# include <time.h>

# define MAX(a, b) (a > b ? a : b)
# define MIN(a, b) (a < b ? a : b)

/*Function*/
double **CREAT_MATRIX(); 
double complex **CREAT_MATRIX_COMPLEX();
void SVL_TRIANGLE ();
void SVL_FFTW ();
void SVL_IFFTW ();
void SVL_SWAPQUADRANTS ();
void SVL_TRUNCATE_FFTW_SPATIAL_HARMONIC ();
//void SVL_GRADING_VECTOR (); 
void SVL_ORIENTATION_FUNCTION ();
void SVL_FDDER();
void Cart2Pol();
void Pol2Cart();
void INVERSE_REAL_MATRIX();
void PRINT_2D_ARRAY_2D();
void PRINT_1D_ARRAY_2D();
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
  int Nx = 128;
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
  
  int m0, n0, n1, n2, m1, m2; 
                          // floor: return the last integral value less than o equal to x
  m0 = 1 + floor(Nx/2);  // x-axis middle point, ex: floor of 12.4 is 12.0
  n0 = 1 + floor(Ny/2);  // y-axis middle point, ex: floor of 12.4 is 12.0

  m1 = m0 - floor(NM/2) - 1; // ex: floor of 12.6 is 12.0
  m2 = m0 + floor(NM/2) - 1; // ex: floor of 13.1 is 13.0
  n1 = n0 - floor(NN/2) - 1; // ex: floor of -2.3 is -3.0
  n2 = n0 + floor(NN/2) - 1; // ex: floor of -3.8 is -4.0
   
  int array[ ] = { m1, m2, n1, n2 };
  int array_size  = 4;

/* Spatial Harmonic - Truncated Array*/
  double complex **TA = CREAT_MATRIX_COMPLEX (NM, NN);
  
  SVL_TRUNCATE_FFTW_SPATIAL_HARMONIC (array, array_size, Nx, Ny, NM, NN, A_re, A_im, TA); 

  printf("\n         Truncated Harmonic Values         \n");
  printf("(m0, m1, m2, n0, n1, n2) = (%d, %d, %d, %d, %d, %d)\n", m0, m1, m2, n0, n1, n2);
  printf("\n         Vector Column AMN(:) \n");
  
/* Spatial Harmonics - Rewritting TA as a Column Array, AMN = TA(:) */
  double complex *AMN  =  malloc(NM * NN * sizeof(double complex)); 
  for (j = 0; j < NN; ++j) {
       for (i = 0; i < NM; ++i) {
	   AMN[i*NM + j]  =   TA[i][j]; 
      }	
   } 
/* Free Memory */
  free(A_re);
  free(A_im); 
  free(TA);
/*************************************************************/
/*                        GRADING VECTOR                     */
/*************************************************************/    
  double pi = 4.0 * atan(1.0);
  double *KX   =  malloc(NM * NN * sizeof(double));
  double *KY   =  malloc(NM * NN * sizeof(double));
    
  printf("\nSTEP 6: (KX, KY) Wavevector cartician coord, writing KX(:) and KY(:) in a column array fashion\n");

/* Writting Grading vectors in Column array fashion, KX(:) and KY(:) */
  for (j = 0; j < NN; ++j) {
      for (i = 0; i < NM; ++i) {
            KX[i*NM + j]  =  (2*pi/Lx)*(j - NM/2); 
	    KY[i*NM + j]  =  (2*pi/Ly)*(i - NN/2);         
        } 
    }
 
  /*SVL_GRADING_VECTOR (pi, NN, NM, Lx, Ly, KX, KY); */
   
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
  int NS[2] = {New_Nx, New_Ny}; 
  int BC[2]  = {1, 1};
  double RES[2] = {New_dx, New_dy};

  int M = New_Nx * New_Ny;          // COMPUTE SIZE OF MATRICES
  
/* Declare the triplet arrays */
   struct cs_sparse *DX  = cs_spalloc(M, M, 2*M, 1, 1);
   struct cs_sparse *D2X = cs_spalloc(M, M, 2*M, 1, 1);
   struct cs_sparse *DY  = cs_spalloc(M, M, 2*M, 1, 1);
   struct cs_sparse *D2Y = cs_spalloc(M, M, 2*M, 1, 1);

/* D = [DX; DY] OPERATOR IS BUILD*/
   struct cs_sparse *D   = cs_spalloc(2*M, M, 4*M, 1, 1);
   SVL_FDDER(M, NS, RES, BC, DX, D2X, DY, D2Y, D);

/* Free memory  */
   cs_spfree(DX);
   cs_spfree(D2X);
   cs_spfree(DY);
   cs_spfree(D2Y);
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
/*									*/
/* cs.entry : it is useful if the number of entries in the matrix is not*/
/*            known when the matrix is first allocated. If space is not */
/*            sufficient for the next entry, the size of the T->i, T->j,*/
/*            and T->x arrays is doubled. The dimensions of T are       */
/*            increased as needed.                                      */
/*									*/
/*   csi cs_entry (cs *T, csi i, csi j, double x)                       */
/*       cs_entry (cs *T, int i, int j, double x)			*/
/* cs_compress : it converts this triplet-form T into a compressed-     */
/*               column matrix C. First, C and a size-n workspace       */ 
/*               are allocated. Next, the number of entries in each     */
/*               column of C is computed, and the column pointer array  */
/*               Cp is constructed as the cumulative sum of the column  */
/*               counts. The counts in w are also replaced with a copy  */
/*              of Cp. cs_compress iterates through each entry in the   */
/*              triplet matrix.                                         */
/*									*/
/*   cs_compress (const cs *T) 						*/
/************************************************************************/

/************************************************************************/  
/* SOLVE LINEAR EQUATION   AX = b  , A  is mxn over determined matrix   */
/*                                   more rows than columns             */ 
/* Least square:  A X = b => A^TA X = A^T b                             */ 
/*                                X = (A^T A)^(-1) (A^T b)               */
/************************************************************************/
/* Convert the Triplet Matrix "D" into Compressed Column Form "A"*/ 
     cs *A, *AT, *C; 
     A = cs_compress(D);
/* Transpose a Sparse matrix AT = A^T  */
     AT = cs_transpose(A, 1);
/* Sparse Matrix Multiplication C = A^T * A */
     C = cs_multiply(AT, A);
/*Free Memory */
     cs_spfree(D);
     cs_spfree(A);
/*************************************************************/
/*                    MAIN LOOP SPATIAL VARIANCE             */
/*************************************************************/
     printf("\nSTEP 9:SPATIAL VARIANCE LOOP\n");
  
/* Output*/  
     double         **PHI = CREAT_MATRIX (M, M); 
     double complex  **UC = CREAT_MATRIX_COMPLEX (M, M);
     double complex   **S = CREAT_MATRIX_COMPLEX (M, M);

     double  *f   = malloc(2 * M * sizeof(double)); 
     double  *B   = malloc(M * sizeof(double));
     double **Kx  = CREAT_MATRIX (New_Nx, New_Ny);
     double **Ky  = CREAT_MATRIX (New_Nx, New_Ny);
     double **TH  = CREAT_MATRIX (New_Nx, New_Ny);
     double **RHO = CREAT_MATRIX (New_Nx, New_Ny);
     int nk;

/*************************************************************/
/*                     START MAIN LOOP                       */
/*************************************************************/
   int NK = NM*NN; // length of AMN

   for (nk = 0; nk < NK; nk++) { 
/****************************************************************/  
/*          BUILDING THE GRADING VECTOR FOR EACH nk             */
/*  K(r) = 2*pi/Delta(r)[a_x*cos(THETA(r)) + a_y*sin(THETA(r))] */
/****************************************************************/  
      //printf(" loop interations nk = %d\n",nk); 
      for (i = 0; i < New_Nx; i++) { 
	  for (j = 0; j < New_Ny; j++) {
            Kx[i][j]  = KX[nk];
	    Ky[i][j]  = KY[nk];  
	  }
       }
/*************************************************************/
/*   Attribute 1 - Lattice Orientation Function  THETA(r)    */
/*              K field is oriented an angle THETA(r)        */ 
/*************************************************************/  

/* TRANSFORMS CARTESIAN COORDINATES TO POLAR COORDINATES, output: TH and RHO */
      Cart2Pol(New_Nx, New_Ny, Kx, Ky, TH, RHO); 

      for (i = 0; i < New_Nx; i++) { 
	for (j = 0; j < New_Ny; j++) {
           TH[i][j] = TH[i][j] + THETA[i][j] ;
	}
      }
 
/* TRANSFORMS THE POLAR COORDINATE TO 2D CARTESIAN COORDINATES, output: Kx and Ky */     
      Pol2Cart(New_Nx, New_Ny, TH, RHO, Kx, Ky); 

/***************************************************************/
/* Attribute 2 - Lattice Spacing Function  RHO = RHO/Delta(r)  */ 
/*              K field is divie by Delta(r)                   */ 
/***************************************************************/
     /*int atri_2 = 1;
     for (i = 0; i < New_Nx; i++) { 
	for (j = 0; j < New_Ny; j++) {
             RHO[i][j] = RHO[i][j]/atri_2 ;
	}
      }*/

/*************************************************************/
/*    Attribute 3 -  Fill Fraction Function  f(r)            */ 
/*************************************************************/
      
/* WRITING f AS COLUMN ARRAY, f = [Kx(:); Ky(:)] each Kx column is order one under */
/* the other, next Ky columns are order on the same way */            
       for (j = 0; j < 2*New_Ny; j++) {
           for (i = 0; i < New_Nx; i++) {                  
	      if (j < New_Ny) {    //New_Nx = First Leading Dimension in Array (LDA), number of rows, A(M,N) then LDA = New_Nx
                f[i + j*New_Ny] = Kx[i][j];	    
	      } else {
	        f[i + j*New_Ny] = Ky[i][j - New_Ny];	    
	      }
          }
       }


/************************************************************************/
/*          Matrix-vector multiplication B = A^T*f + B,                */
/************************************************************************/
       cs_gaxpy (AT, f, B); // return B, carefull B will add , make sure it will be update to zeros
     
/************************************************************************/
/* Least square  A X = b => (A^T A) X = (A^T b) => C X = B              */ 
/************************************************************************/
/************************************************************************/
/* csi cs_lusol (csi order, const cs *A, double *b, double tol)         */
/* X = C\B where C is unsymmetric; B overwritten with solution          */
/* Note that the vector B overwritten to contain the solution vector X  */ 
/* (i.e return B)                                                       */
/************************************************************************/
     cs_lusol(1, C, B, 1); // return B time = 19 s
     //cs_qrsol(3, C, B); // return B  time = 91 s
     //cs_cholsol(0, C, B); // return B  time = 51 s

/*Rewrite solution in a 2D array of shape (Nx_New x Ny_Newy) */
     for (j = 0; j < New_Ny; j++) {
         for (i = 0; i < New_Nx; i++) {      
             PHI[i][j] = B[i + New_Nx*j]; 
             B[i + New_Nx*j] = 0; //clean B, 
        }
     }
     
/*************************************************************/    
/*           SOLVE the exponential S = AMN[nk]*exp(i*PHI)    */
/*************************************************************/
     for (i = 0; i < New_Nx; i++) { 
          for (j = 0; j < New_Ny; j++) {                      
	      S[i][j]   = AMN[nk] * cexp(I*PHI[i][j]);
              UC[i][j]  = UC[i][j] + S[i][j]; 
          }
      }   
      
     if (nk == 10) {
        FILE *fp11, *fp12, *fp13; // open a file
        fp11 = fopen("data11", "w"); // PHI
        fp12 = fopen("data12", "w"); // S
        fp13 = fopen("data13", "w"); // UC

        for (i = 0; i < New_Nx; i++) {
            for (j = 0; j < New_Ny; j++) {
                fprintf(fp11,"%f   ",PHI[i][j]); 
	        fprintf(fp12,"%f   ",creal(S[i][j]));  
                fprintf(fp13,"%f   ",creal(UC[i][j]));  
            }  
            fprintf(fp11, "\n");
            fprintf(fp12, "\n");	
            fprintf(fp13, "\n");		
        }
        fclose(fp11);
        fclose(fp12);
        fclose(fp13);
      }
 } /* end nk array */


/*Clean up memory usage.*/ 
  cs_spfree(AT);
  cs_spfree(C);

  free(f);
  free(B);
  free(Kx);
  free(Ky);
  free(TH);
  free(RHO);

  FILE *fp14, *fp15, *fp16; // open a file
  fp14 = fopen("data14", "w"); // PHI
  fp15 = fopen("data15", "w"); // S
  fp16 = fopen("data16", "w"); // UC
  
  for(i = 0; i < New_Nx; i++) {
      for(j = 0; j < New_Ny; j++) {
          fprintf(fp14,"%f   ",PHI[i][j]); 
	  fprintf(fp15,"%f   ",creal(S[i][j]));  
          fprintf(fp16,"%f   ",creal(UC[i][j]));  
    }  
    fprintf(fp14, "\n");
    fprintf(fp15, "\n");	
    fprintf(fp16, "\n");		
   }

  fclose(fp14);
  fclose(fp15);
  fclose(fp16);

   clock_t toc = clock(); /*End Elapsed time */
   printf("\n   Elapsed: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);

  printf("\n=============THE END===========\n");
  return 0;
}
/*************************************************************/
/*                      End main program                     */
/*************************************************************/
/*************************************************************/
/*      INPUT 1D ARRAY (==> PRINT ==>) OUTPUT 2D ARRAY       */
/*************************************************************/
void PRINT_1D_ARRAY_2D( char* desc, int m, int n, double* a, int lda ) {
        int i, j; //LDA = first dimension in array, number of rows, A(M,N) then LDA=M
        printf( "\n %s\n", desc);
        for( i = 0; i < m; i++ ) {
            for( j = 0; j < n; j++ ) {
               printf("%2.1f  ",a[i+j*lda] );
            }
            printf( "\n" );
        }
}/* END PRINT_1D_ARRAY_2D  */  
/*************************************************************/
/*     INPUT 2D ARRAY (==> PRINT ==>) OUTPUT 2D ARRAY        */
/*************************************************************/
void PRINT_2D_ARRAY_2D(char* desc, int m, int n, double **a) {    
   int i, j;

   printf( "\n %s\n", desc );

   for (i = 0; i < m; ++i) {
       for (j = 0; j < n; ++j) { 
           printf("%2.1f  ",a[i][j]);
       }
      printf("\n");	
   }
} /* END PRINT_2D_ARRAY_2D  */  
/*************************************************************/
/*                   CARTESIAN TO POLAR                     */
/*************************************************************/ 
void Cart2Pol(int New_Nx, int New_Ny, double **Kx, double **Ky, double **TH ,double **RHO){ 

  int i, j;

  for (i = 0; i < New_Nx; ++i) {
     for (j = 0; j < New_Ny; ++j)  {
         RHO[i][j] = sqrt(Kx[i][j] * Kx[i][j] + Ky[i][j] * Ky[i][j] );
         TH[i][j] =  atan2(Ky[i][j], Kx[i][j]);
     }    
  }
 }  /*   END CART2POL     */
/*************************************************************/
/*                 POLAR TO CATRTESIAN                       */
/*************************************************************/
void Pol2Cart(int New_Nx, int New_Ny, double **TH, double **RHO, double **Kx, double **Ky) { 

  int i, j;
 
  for (i = 0; i < New_Nx; ++i) {
     for (j = 0; j < New_Ny; ++j)  {
         Kx[i][j] = RHO[i][j] * cos(TH[i][j]);
         Ky[i][j] = RHO[i][j] * sin(TH[i][j]);
     }    
  }
}  /*     END POL2CART     */
/*************************************************************/
/*                 Call Creat Matrix                       */
/*************************************************************/
 double **CREAT_MATRIX(int rows, int cols) {
 int j; 
 double **mat; /* Defining a temporal "mat" array , otherwise would have to use (*memory) everywhere U is used (yuck) */

/* Each row should only contain double*, not double**, because each row will be an array of double */
 mat = (double **) malloc(rows * sizeof(double *)); // create Nx-row temporal pointers array 

 if (mat == NULL) {
    printf("Failure to allocate room for row pointers.\n ");
    exit(0);
 }

/* Had an error here.  Alloced rows above so iterate through rows not cols here */
   for (j = 0; j < rows; j++) {
/* Allocate array, store pointer  */
     mat[j] = (double *) malloc(cols * sizeof(double)); // Create on each temporal rows we allocate Ny-columns temporal pointer
    
      if (mat[j] == NULL) {
            printf("Failure to allocate for row[%d]\n", j);
            exit(0); 
/* still a problem here, if exiting with error, should free any column mallocs that were successful. */
       } 
 } 

  return (mat); // return mat as U
  free(mat); // free mat No U
} /* END FUNCTION CREAT MATRIX */
/*************************************************************/
/*                 Call Creat Matrix Complex                        */
/*************************************************************/
 double complex **CREAT_MATRIX_COMPLEX(int rows, int cols) {
 int j; 
 double complex **mat; /* Defining a temporal "mat" array , otherwise would have to use (*memory) everywhere U is used (yuck) */

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

/*************************************************************/
/*                  ADD TRIANGLE                             */
/* INSERTING AN EQUILATERAL TRIANGLE INTO A SQUARE UNIT CELL */
/*************************************************************/
void SVL_TRIANGLE (int Nx, int Ny, double **U){
 int nx1, nx2, ny1, ny2;
 int nx, ny, nxa, nxb;
 int i, j;
 double f;
 FILE *fp2; // open a file
 fp2 = fopen("data2", "w"); // save triangle into a Square Unit Cell 
 printf("\nSTEP 1: FILL U WITH ONES TO BUILD TRIANGLE-DEVICE\n");
 
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

/***************************************************/
/*                  FFTW                           */
/*     FAST FOURIER TRANSFORM ON THE WEST          */
/***************************************************/
  void SVL_FFTW(int Nx, int Ny, double **U, double **A_re, double **A_im) {
  int i, j;
   
  FILE *fp3, *fp4; // open a file
  fp3 = fopen("data3", "w"); // save Real
  fp4 = fopen("data4", "w"); // save Imagynary
  printf("\nSTEP 2: SET FFTW DATA\n\n");
  
  fftw_plan My_plan_fft;
// pointers to input and output data 
  fftw_complex  *U_in , *U_out; 
  
// Allocate input arrays
   U_in  = fftw_malloc(sizeof(fftw_complex) * Nx * Ny);
   
// Allocate output arrays
   U_out = fftw_malloc(sizeof(fftw_complex) * Nx * Ny);

// input data fft
 for (i = 0; i < Nx; ++i) {
     for (j = 0; j < Ny; ++j) {
      U_in[i*Ny + j][0] = U[i][j]; // real part
      U_in[i*Ny + j][1] = 0.0;     // imaginary part
      } 
 }

// Create plans for fftw 2D. Here, we set the type of transformation we want, i.e. 1D, 2D or 3D and fft(FORWARD) or ifft(BACKWARD)
   My_plan_fft = fftw_plan_dft_2d(Nx, Ny, U_in, U_out, FFTW_FORWARD, FFTW_ESTIMATE);

// Execute fftw plan, output file U_out 
   fftw_execute(My_plan_fft);

// Normalize values 
   double normalization = (double) (Nx * Ny); 
   
   for (i = 0; i < Nx * Ny; i++) {
       U_out[i][0] = U_out[i][0] / normalization;
       U_out[i][1] = U_out[i][1] / normalization; 
   }

// output data fft  
   for (i = 0; i < Nx; ++i) {
       for (j = 0; j < Ny; ++j) {
           A_re[i][j] =  U_out[i*Ny + j][0];
           A_im[i][j] =  U_out[i*Ny + j][1];
       }
   }

// Destroy plan 
  fftw_destroy_plan(My_plan_fft);  
  fftw_free(U_in);		   //Free memory
  fftw_free(U_out);		   //Free memory
  
  printf("      SVL_FFTW - A_re array is pass by reference from SVL_FFTW\n"); 
  for(i = 0; i < Nx; i++) {
    for(j = 0; j < Ny; j++) {
       fprintf(fp3,"%f   ",A_re[i][j]);
       //printf("%f   ",A_re[i][j]);
    }
    fprintf(fp3,"\n");
    //printf("\n");
  }
  
  printf("\n      SVL_FFTW - A_im array is pass by reference from SVL_FFTW\n");
  for(i = 0; i < Nx; i++) {
    for(j = 0; j < Ny; j++) {
       fprintf(fp4,"%f   ",A_im[i][j]);
       //printf("%f   ",A_im[i][j]);
    }    
    fprintf(fp4,"\n");
    //printf("\n");
  } 

  fclose(fp3);
  fclose(fp4);

} /* END FUNCTION SL_FFTW */

/***************************************************/
/*                  IFFTW                           */
/***************************************************/
void SVL_IFFTW (int Nx, int Ny, double **A_re, double **A_im, double **inv_A_re, double **inv_A_im) {

  int i, j;
  printf("\nSTEP 3: SET IFFTW DATA\n"); 
    
  fftw_plan My_plan_ifft;
  fftw_complex *inv_U_in, *inv_U_out;  // pointers to input and output data
 
// Allocate input arrays
  inv_U_in  = fftw_malloc(sizeof(fftw_complex) * Nx * Ny);
  
// Allocate output arrays
  inv_U_out = fftw_malloc(sizeof(fftw_complex) * Nx * Ny);

// Creat ifftw plan for fftw 2D. Here, we set what kind of transformation we want, i.e. 1D, 2D or 3D and fft(FORWARD) or ifft(BACKWARD)
  My_plan_ifft = fftw_plan_dft_2d(Nx, Ny, inv_U_in, inv_U_out, FFTW_BACKWARD, FFTW_ESTIMATE);

  for (i = 0; i < Nx; ++i) {
      for (j = 0; j < Ny; ++j) {
          inv_U_in[i*Ny + j][0] = A_re[i][j]; // real part
          inv_U_in[i*Ny + j][1] = A_im[i][j]; // imaginary part
      } 
  }

// Execute ifftw plan, output file inv_U_out 
   fftw_execute(My_plan_ifft);

//output data ifft
/*   for (i = 0; i < Nx; ++i) {
       for (j = 0; j < Ny; ++j) {
           inv_A_re[i][j] =  inv_U_out[i*Ny + j][0];
           inv_A_im[i][j] =  inv_U_out[i*Ny + j][1];
      //printf("{ %g, %g }  ", inv_A_re[i][j], inv_A_im[i][j]);
       }
    //printf("\n");
   }
*/
// Destroy  inv plan 
  fftw_destroy_plan(My_plan_ifft);  
  fftw_free(inv_U_in);		   //Free memory
  fftw_free(inv_U_out);		   //Free memory
} /* END FUNCTION SL_IFFTW */

/***************************************************/
/* SWAP QUADRANTS (1 <--> 3, 2 <--> 4) DIAGONALLY  */
/***************************************************/
void SVL_SWAPQUADRANTS (int Nx, int Ny, double **A_re, double **A_im) {

  int i, j;
  int N2y, N2x; 
  double tmp13, tmp24;
   
  FILE *fp5, *fp6; // open a file
  fp5 = fopen("data5", "w"); // save Real
  fp6 = fopen("data6", "w"); // save Imagynary
  printf("\nSTEP 4: SET FFTW_SWAP\n");
  
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

/***************************************************/
/*     TRUCATE THE SPATIAL HARMONICS               */
/***************************************************/
 void SVL_TRUNCATE_FFTW_SPATIAL_HARMONIC (int array[ ], int arraySize, int Nx, int Ny,
	      int NM, int NN, double **A_re, double **A_im, double complex **TA) { 
// notice the function is void as it returns nothing
// you passed in the address of the variables
// pointers are a pointer to an address
// these statements are basically saying that the address contains the following value
  int i, j;
  printf("\nSTEP 5: SET TRUNCATION ARRAY OF FFTW SPATIAL HARMONICS\n");
   
  FILE *fp7, *fp8; //open a file
  fp7 = fopen("data7", "w"); 
  fp8 = fopen("data8", "w");  
  int n1, n2, m1, m2;  

  m1 = array[0]; // ex: floor of 12.6 is 12.0
  m2 = array[1]; // ex: floor of 13.1 is 13.0
  n1 = array[2]; // ex: floor of -2.3 is -3.0
  n2 = array[3]; // ex: floor of -3.8 is -4.0
  
  
// MAGNITUDE
  for (i = m1 ; i <= m2; ++i) {
      for (j = n1 ; j <= n2; ++j) {
	   //TA[i-m1][j-n1] = sqrt((A_re[i][j] * A_re[i][j]) + (A_im[i][j] * A_im[i][j]));
           TA[i-m1][j-n1] = A_re[i][j] + I*A_im[i][j];
           fprintf(fp7,"%f   ",  creal(TA[i-m1][j-n1])); 
           fprintf(fp8,"%f   ",  cimag(TA[i-m1][j-n1]));           
       }	
       fprintf(fp7,"\n");
       fprintf(fp8,"\n");
   } 

   fclose(fp7); // close a file
   fclose(fp8); // close a file
} //END FUNCTION TRUNCATE FFTW SPATIAL HARMONIC 

/***************************************************/
/*          SVL_GRADING_VECTOR                     */
/***************************************************/
/*
 This part is not working 
 void SVL_GRADING_VECTOR (double pi, int NN, int NM, int Lx, int Ly, int *KX, int *KY)  {
 */

/***************************************************/
/*          SVL ORIENTATION FUNCTION               */
/***************************************************/
void SVL_ORIENTATION_FUNCTION (double pi, int NN, int NM, int Lx, int Ly, int NPx, int NPy,
        int New_Nx, int New_Ny, double New_dx, double New_dy, double **RSQ, double **THETA) {
  
   
  // Calculate grating vectores 
  int i, j;
  double xa, ya, mean_xa, mean_ya;
  double X[New_Nx][New_Ny], Y[New_Nx][New_Ny];
  
  printf("\nSTEP 7: RSO and  THETA\n");
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

/***************************************************
                    SVL FDDER                    
 FDDDER      Finite-Difference Derivative Operators
 *************************************************** 
 [DX,D2X,DY,D2Y] = fdder(NS,RES,BC);

 This function generates matrix derivative operators
 for scalar functions on a collocated grid.

 INPUT ARGUMENTS
 ================
 NS        [Nx Ny] size of grid
 RES       [dx dy] grid resolution
 BC        [BCx BCy] Boundary Conditions
           0 = Dirichlet, -1 = Periodic, +1 = Neumann

 OUTPUT ARGUMENTS
 ================
 DX        First-order derivative with respect to x
 D2X       Second-order derivative with respect to x
 DY        First-order derivative with respect to y
 D2Y       Second-order derivative with respect to y

 ************************************************************************/
/*   D[m][n] is a 2D operator of dimension (2*M x M)                    */ 
/*   nzmax : (int) maximum number of entries 				*/
/*   m  : (int) number of rows 	of the operator			        */
/*   n  : (int) number of columns of the operator			*/
/*   *p : (int) column pointers (size n+1) or col indices (size nzmax) 	*/	
/*   *i : (int) row indices, size nzmax 				*/	
/*   *j : (int) 							*/
/*   *x : (double) numerical values, size nzmax 			*/
/*   nz : (int) # of entries in triplet matrix, -1 for compressed-col	*/
/*									*/
/* allocate a sparse matrix (triplet form or compressed-column form)    */
/*   cs *cs_spalloc (csi m, csi n, csi nzmax, csi values, csi triplet)  */
/*   cs_spalloc (int m, int n, int nzmax, int values, int triplet)	*/
/*									*/
/* add an entry to a triplet matrix; return 1 if ok, 0 otherwise        */
/*   csi cs_entry (cs *T, csi i, csi j, double x)                       */
/*   cs_entry (cs *T, int i, int j, double x)				*/
/*   cs_compress (const cs *T) 						*/
/************************************************************************/

void SVL_FDDER(int M, int NS[2], double RES[2], int BC[2], 
              struct cs_sparse *DX, struct cs_sparse *D2X, struct cs_sparse *DY, struct cs_sparse *D2Y, struct cs_sparse *D) {
  int i, k;
  double CX, C2X, CY, C2Y;

// DETERMINE SIZE OF GRID
  int Nx = NS[0];
  int Ny = NS[1];

  printf("\nSTEP 8: SET FDDER\n");
  printf(" BUILD THE STANDART-TRIPLET - THREE COLUMN ARRAY");

/* counter */
int sx1 = 0, sx2 = 0;

/* Unit cells Arrays */
double **DXC, **D2XC, *DY1, *DY2; 

DXC  = CREAT_MATRIX (Nx * Nx, 3);
D2XC = CREAT_MATRIX (Nx * Nx, 3);

DY1  = (double *) malloc(Nx * sizeof(double));
DY2  = (double *) malloc(Nx * sizeof(double));

/******************************************************/
/*                   DX and D2X                       */
/* CREATE MATRIX IF DIMENSION EXISTS (True if Nx > 0) */
/******************************************************/
/* Build Unit cells array  DXC and D2XC*/
if (Nx > 0 ) { 
// Build Default DX D2X
   for (i = 0; i < Nx - 2; i++) {
       //DXC[i + 1][i] = 1; 
       DXC[sx1][0] = i + 1;
       DXC[sx1][1] = i;
       DXC[sx1][2] = -1; 
       //DXC[i + 1][i + 2] =-1;
       DXC[sx1 + 1][0] = i + 1;
       DXC[sx1 + 1][1] = i + 2;
       DXC[sx1 + 1][2] = 1;
       sx1 = sx1 + 2; 
   }

   for (i = 0; i < Nx - 2; i++) {
       //D2XC[i + 1][i] = 1; 
       D2XC[sx2][0] = i + 1;
       D2XC[sx2][1] = i;
       D2XC[sx2][2] = 1;
       //D2XC[i + 1][i + 2] = 1;
       D2XC[sx2 + 1][0] = i + 1;
       D2XC[sx2 + 1][1] = i + 2;
       D2XC[sx2 + 1][2] = 1;
       sx2 = sx2 + 2;
   }
 
   for (i = 1; i < Nx - 1 ; i++) {
       //D2XC[i][i] = -2;
       D2XC[sx2][0] = i;
       D2XC[sx2][1] = i;
       D2XC[sx2][2] = -2;
       sx2++;
   }

// Incorporate Boundary Conditions
    switch (BC[0]) {
        case -1: //Periodic
                //DXC[0][Ny - 1]  = -1;
                DXC[sx1][0] = 0;
                DXC[sx1][1] = Ny - 1;
                DXC[sx1][2] = 1;
                //DXC[Nx - 1][0]  =  1;
       		DXC[sx1 + 1][0] = Nx - 1;
       		DXC[sx1 + 1][1] = 0;
       		DXC[sx1 + 1][2] = 1;
                sx1 = sx1 + 2;
                //D2XC[0][Ny - 1] =  1;
                D2XC[sx2][0] = 0;
                D2XC[sx2][1] = Ny - 1;
                D2XC[sx2][2] = 1;
                //D2XC[Nx - 1][0] =  1;
                D2XC[sx2 + 1][0] = Nx - 1;
                D2XC[sx2 + 1][1] = 0;
                D2XC[sx2 + 1][2] = 1;
                sx2 = sx2 + 2;
    
        case +1: //Neumann
                //DXC[0][0] = -2;
                DXC[sx1][0] = 0;
                DXC[sx1][1] = 0;
                DXC[sx1][2] = -2;
                //DXC[0][1] =  2;
                DXC[sx1 + 1][0] = 0;
                DXC[sx1 + 1][1] = 1;
                DXC[sx1 + 1][2] = 2;
                //DXC[Nx - 1][Ny - 1] =  2;
                DXC[sx1 + 2][0] = Nx - 1;
                DXC[sx1 + 2][1] = Ny - 2;
                DXC[sx1 + 2][2] = -2;
                //DXC[Nx - 1][Ny - 2] = -2;
                DXC[sx1 + 3][0] = Nx - 1;
                DXC[sx1 + 3][1] = Ny - 1;
                DXC[sx1 + 3][2] = 2;
                sx1 = sx1 + 4;

       } // end swicth

// Finish Calculation
    CX  = (double) 1/(2*RES[0]);
    C2X = (double) 1/(RES[0]*RES[0]);

    /*for (k = 0; k < Nx; k++) {
        for (i = 0; i < sx1; i++) {
             if (DXC[i][2] != 0) { cs_entry(DX, DXC[i][0] + k*Nx, DXC[i][1] + k*Nx, CX*DXC[i][2]);}
    	}
    }*/
    
    for (k = 0; k < Nx; k++) {
        for (i = 0; i < sx2; i++) {
            if (D2XC[i][2] != 0) { cs_entry(D2X, D2XC[i][0] + k*Nx, D2XC[i][1] + k*Nx, C2X*D2XC[i][2]); }
    	}
    }

/* Main operator D = [DX ; DY] */
    for (k = 0; k < Nx; k++) {
        for (i = 0; i < sx1; i++) {
             if (DXC[i][2] != 0) { cs_entry(D, DXC[i][0] + k*Nx, DXC[i][1] + k*Nx, CX*DXC[i][2]);}
    	}
    }
 }  // end if Nx   
/******************************************************/
/*                     DY and D2Y                     */
/* CREATE MATRIX IF DIMENSION EXISTS (True if Ny > 0) */
/******************************************************/
/* Build Unit cells array  DY1 and DY2*/
 if (Ny > 0) {
   for (i = 0; i < Nx; i++) { 
       DY1[i]  = 1;           
       DY2[i]  = 2;
   }
   
   CY  = (double) 1/(2*RES[1]);
   C2Y = (double) 1/(RES[1]*RES[1]);

// Incorporate Boundary Conditions
  switch (BC[1]){
        case -1: //periodic
             for (i = 0; i < Nx; i++) {
                 //DY[i + (Nx - 1)*Nx][i] = DY1[i][i];
                 //DY[i][i + (Nx - 1)*Nx] = DY1[i][i];
                 //cs_entry(DY, i + (Nx - 1)*Nx, i, CY * DY1[i]);
                 //cs_entry(DY, i, i + (Nx - 1)*Nx,  CY * DY1[i]);

                 /* Main operator D = [DX ; DY] */
                 cs_entry(D, M + i + (Nx - 1)*Nx, i, CY * DY1[i]);
                 cs_entry(D, M + i, i + (Nx - 1)*Nx,  CY * DY1[i]);

                 //D2Y[i + (Nx - 1)*Nx][i] = DY1[i][i];
                 //D2Y[i][i + (Nx - 1)*Nx] = DY1[i][i];
                 cs_entry(D2Y, i + (Nx - 1)*Nx, i,  C2Y * DY1[i]);
                 cs_entry(D2Y, i, i + (Nx - 1)*Nx,  C2Y * DY1[i]);
              }
        case +1: //extended
              for (i = 0; i < Nx; i++) {
                  //DY[i][i] =-DY2[i][i];
                  //DY[i][i + Ny] = DY2[i][i];
                  //cs_entry(DY, i , i, - CY * DY2[i]);
                  //cs_entry(DY, i, i + Ny,  CY * DY2[i]);

                  /* Main operator D = [DX ; DY] */
                  cs_entry(D, M + i , i, -CY * DY2[i]);
                  cs_entry(D, M + i, i + Ny,  CY * DY2[i]);
               }

               for (k = 1; k < Nx - 1; k++) {
                   for (i = 0; i < Nx; i++) {
                      //DY[i + k*Nx][i + (k - 1)*Nx] = DY1[i][i];
                      //DY[i + k*Nx][i + (k + 1)*Nx] = DY1[i][i];
                      //cs_entry(DY, i + k*Nx, i + (k - 1)*Nx, - CY * DY1[i]);
                      //cs_entry(DY, i + k*Nx, i + (k + 1)*Nx,  CY * DY1[i]);

                      /* Main operator D = [DX ; DY] */
                      cs_entry(D, M + i + k*Nx, i + (k - 1)*Nx, -CY * DY1[i]);
                      cs_entry(D, M + i + k*Nx, i + (k + 1)*Nx,  CY * DY1[i]);
    	           }
               }

               for (k = 1; k < Nx - 1; k++) {
                   for (i = 0; i < Nx; i++) {
                      //D2Y[i + k*Nx][i + (k - 1)*Nx] = DY1[i][i];
                      //D2Y[i + k*Nx][i + k*Nx] = -DY2[i][i];
                      //D2Y[i + k*Nx][i + (k + 1)*Nx] = DY1[i][i];                     
                      cs_entry(D2Y, i + k*Nx, i + (k - 1)*Nx,  C2Y * DY1[i]);
                      cs_entry(D2Y, i + k*Nx, i+ k*Nx,- C2Y * DY2[i]);
                      cs_entry(D2Y, i + k*Nx, i + (k + 1)*Nx,  C2Y * DY1[i]);
    	           }
               }

               for (i = 0; i < Nx; i++) {
                   //DY[i + (Nx - 1)*Nx][i + (Nx - 2)*Nx] =-DY2[i][i];
                   //DY[i + (Nx - 1)*Nx][i + (Nx - 1)*Nx] = DY2[i][i];
                   //cs_entry(DY, i + (Nx - 1)*Nx, i + (Nx - 2)*Nx,- CY * DY2[i]);
                   //cs_entry(DY, i + (Nx - 1)*Nx, i + (Nx - 1)*Nx,  CY * DY2[i]);

                   /* Main operator D = [DX ; DY] */
                   cs_entry(D, M + i + (Nx - 1)*Nx, i + (Nx - 2)*Nx, -CY * DY2[i]);
                   cs_entry(D, M + i + (Nx - 1)*Nx, i + (Nx - 1)*Nx,  CY * DY2[i]);
               }
    } // end switch

 } //end if Ny 

 } /* END SVL FDDER  */




