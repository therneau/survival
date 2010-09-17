/* $Id: chinv2.c 11357 2009-09-04 15:22:46Z therneau $
**
** matrix inversion, given the FDF' cholesky decomposition
**
** input  **matrix, which contains the chol decomp of an n by n
**   matrix in its lower triangle.
**
** returned: the upper triangle + diagonal contain (FDF')^{-1}
**            below the diagonal will be F inverse
**
**  Terry Therneau
*/
#include "survS.h"
#include "survproto.h"

void chinv2(double **matrix , int n)
     {
     register double temp;
     register int i,j,k;

     /*
     ** invert the cholesky in the lower triangle
     **   take full advantage of the cholesky's diagonal of 1's
     */
     for (i=0; i<n; i++){
	  if (matrix[i][i] >0) {
	      matrix[i][i] = 1/matrix[i][i];   /*this line inverts D */
	      for (j= (i+1); j<n; j++) {
		   matrix[j][i] = -matrix[j][i];
		   for (k=0; k<i; k++)     /*sweep operator */
			matrix[j][k] += matrix[j][i]*matrix[i][k];
		   }
	      }
	  }

     /*
     ** lower triangle now contains inverse of cholesky
     ** calculate F'DF (inverse of cholesky decomp process) to get inverse
     **   of original matrix
     */
     for (i=0; i<n; i++) {
	  if (matrix[i][i]==0) {  /* singular row */
		for (j=0; j<i; j++) matrix[j][i]=0;
		for (j=i; j<n; j++) matrix[i][j]=0;
		}
	  else {
	      for (j=(i+1); j<n; j++) {
		   temp = matrix[j][i]*matrix[j][j];
		   if (j!=i) matrix[i][j] = temp;
		   for (k=i; k<j; k++)
			matrix[i][k] += temp*matrix[j][k];
		   }
	      }
	  }
     }
