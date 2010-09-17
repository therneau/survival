/* $Id: chinv3.c 11166 2008-11-24 22:10:34Z therneau $ */
/*
** matrix inversion, given the cholesky decomposition
** This is a specialized form for the frailty problem.  The matric C in this
**   case has C[1:m, 1:m] diagonal and  C[(m+1):n, 1:n)] is dense. 
**
** arguments are:
**     n         the size of the matrix to be factored
**     m         the size of the diagonal upper portion
**     diag      the diagonal upper portion
**     **matrix  a ragged array containing the dense portion
**
** returned: the lower triange will contain the inverse of the cholesky
**
**  Terry Therneau
*/
#include "survS.h"
#include "survproto.h"

void chinv3(double **matrix , int n, int m, double *fdiag)
     {
     int i,j,k;
     int n2, ii;
   
     n2 = n-m;   /* number of full covariates */

     /*
     ** invert the cholesky in the lower triangle
     **   take full advantage of the cholesky's diagonal of 1's
     ** start with the diagonal part
     */
     for (i=0; i<m; i++) {
	 if (fdiag[i] >0){
	     fdiag[i] = 1/fdiag[i];     /* this line inverts D */
	     for (j=0; j<n2; j++) 
		 matrix[j][i] = -matrix[j][i];
	     }
	  }
	 
     /* Now for the original portion */
     for (i=0; i<n2; i++){
	  ii = i+m;
	  if (matrix[i][ii] >0) {
	      matrix[i][ii] = 1/matrix[i][ii];   /*this line inverts D */
	      for (j= (i+1); j<n2; j++) {
		   matrix[j][ii] = -matrix[j][ii];
		   for (k=0; k<ii; k++)     /*sweep operator */
			matrix[j][k] += matrix[j][ii]*matrix[i][k];
		   }
	      }
	  }
     }

/*
** Do a specialized matrix product
**  This is split out of chinv3 so that one more number may be saved by
**  the calling routine.
*/
void chprod3(matrix, n, m, fdiag)
int  n, m;
double fdiag[];
double **matrix;
     {
     double temp;	 
     int i,j,k;
     int n2, ii;
   
     n2 = n-m;   /* number of full covariates */
     /*
     ** lower triangle now contains inverse of cholesky
     ** calculate F'DF (inverse of cholesky decomp process) to get inverse
     **   of original matrix
     **make chtest3 

     ** Although the original matrix is sparse, and so are the cholesky
     **  factors, the inverse of the original is not sparse!
     **  But we want to multiply out only the corner anyway
     */
     for (i=0; i<n2; i++) {
	  ii = i+m;
	  if (matrix[i][ii]==0) {  /* singular row */
		for (j=0; j<i; j++) matrix[j][ii]=0;
		for (j=ii; j<n; j++) matrix[i][j]=0;
		}
	  else {
	      for (j=(i+1); j<n2; j++) {
		   temp = matrix[j][ii]*matrix[j][j+m];
		   if (j!=i) matrix[i][j+m] = temp;
		   for (k=i; k<j; k++)
			matrix[i][k+m] += temp*matrix[j][k+m];
		   }
	      }
	  }
     }

