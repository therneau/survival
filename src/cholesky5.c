/*
** subroutine to do a generalized Cholesky decompostion on a matrix: C = FDF'
**   where F is lower triangular with 1's on the diagonal, and D is diagonal
** If D is all >0, then C was symmetric positive definite, if D >=0, C is
**  non-negative definite.
**
** The only difference between this routine and cholesky2 is what it does
**  with negative pivots: cholesky2 considers them to be zero.
**
** arguments are:
**     n         the size of the matrix to be factored
**     **matrix  a ragged array containing an n by n submatrix to be factored
**     toler     the threshold value for detecting "singularity"
**
**  The factorization is returned in the lower triangle, D occupies the
**    diagonal and the upper triangle is left undisturbed.
**
**  Return value:  the rank of the matrix 
**
**  If a column is deemed to be redundant, then that diagonal is set to zero.
**
**   Terry Therneau
*/
#include <math.h>
int cholesky5(double **matrix, int n, double toler)
    {
    double temp;
    int  i,j,k;
    double eps, pivot;
    int rank;

    eps =0;
    for (i=0; i<n; i++) {
	if (fabs(matrix[i][i]) > eps)  eps = fabs(matrix[i][i]);
	}
    if (eps==0) eps = toler;
    else eps *= toler;

    rank =0;
    for (i=0; i<n; i++) {
	pivot = matrix[i][i];
	if (isfinite(pivot)== 0 || fabs(pivot) < eps) {
	    for (j=i; j<n; j++) matrix[j][i] =0;  /* zero the column */
	    }
	else  {
	    rank++;
	    for (j=(i+1); j<n; j++) {
		temp = matrix[j][i]/pivot;
		matrix[j][i] = temp;
		matrix[j][j] -= temp*temp*pivot;
		for (k=(j+1); k<n; k++) matrix[k][j] -= temp*matrix[k][i];
		}
	    }
	}
    return(rank);
    }
/*
** matrix inversion, given the FDF' cholesky decomposition
**  The only difference between this and chinv2 is that it allows negative
**  elements in D, and forces zeros in for redundant columns.
**
** input  **matrix, which contains the chol decomp of an n by n
**   matrix in its lower triangle.
**
** returned: the upper triangle + diagonal contain (FDF')^{-1}
**            below the diagonal will be F inverse
**
**  Terry Therneau
*/
void chinv5(double **matrix , int n, int flag)
     {
     double temp;
     int i,j,k;

     /*
     ** invert the cholesky in the lower triangle
     **   take full advantage of the cholesky's diagonal of 1's
     */
     for (i=0; i<n; i++){
	  if (matrix[i][i] != 0) {
	      matrix[i][i] = 1/matrix[i][i];   /*this line inverts D */
	      for (j= (i+1); j<n; j++) {
		   matrix[j][i] = -matrix[j][i];
		   for (k=0; k<i; k++)     /*sweep operator */
			matrix[j][k] += matrix[j][i]*matrix[i][k];
		   }
	      }
	  else {
	      for (j= (i+1); j<n; j++) matrix[j][i] =0;
	      }
	  }

     if (flag==1) return;  
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
