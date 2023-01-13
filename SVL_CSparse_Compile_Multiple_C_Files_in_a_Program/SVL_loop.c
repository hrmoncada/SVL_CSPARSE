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

/*************************************************************/
/*                     START MAIN LOOP                       */
/*************************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <complex.h>
# include "cs.h"
# include "SVL_header.h"

void SVL_LOOP (int M, int New_Nx, int New_Ny, int NM, int NN, double complex *AMN, double *KX, double *KY, double **THETA, struct cs_sparse *D){ 

   int i, j, nk;
   int NK = NM*NN; // length of AMN
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

   printf("\nSTEP 9:SPATIAL VARIANCE LOOP\n");
   //printf("%s:%s:%d \n", __FILE__, __func__, __LINE__);

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
}
