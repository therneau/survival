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
	if (fabs(pivot) < eps) {
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
