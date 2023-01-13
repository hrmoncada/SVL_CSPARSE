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
# include <stdio.h>
# include <stdlib.h>
//# include <complex.h>
# include "SVL_header.h"
# include "cs.h"


void SVL_FDDER(int M, int NS[2], double RES[2], int BC[2],  struct cs_sparse *D) {
  int i, k;
  double CX, C2X, CY, C2Y;

// DETERMINE SIZE OF GRID
  int Nx = NS[0];
  int Ny = NS[1];

  
/* Declare the triplet arrays */
   struct cs_sparse *DX  = cs_spalloc(M, M, 2*M, 1, 1);
   struct cs_sparse *D2X = cs_spalloc(M, M, 2*M, 1, 1);
   struct cs_sparse *DY  = cs_spalloc(M, M, 2*M, 1, 1);
   struct cs_sparse *D2Y = cs_spalloc(M, M, 2*M, 1, 1);

  printf("\nSTEP 8: SET FDDER - BUILD THE STANDART-TRIPLET - THREE COLUMN ARRAY\n");
//printf("%s:%s:%d \n", __FILE__, __func__, __LINE__);

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

/* Free memory  */
   cs_spfree(DX);
   cs_spfree(D2X);
   cs_spfree(DY);
   cs_spfree(D2Y);
 } /* END SVL FDDER  */

