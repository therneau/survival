/* $Id: cholesky3.c 11166 2008-11-24 22:10:34Z therneau $ */
/*
** subroutine to do Cholesky decompostion on a matrix: C = FDF'
**   where F is lower triangular with 1's on the diagonal, and D is diagonal
** This is a specialized form for the frailty problem.  The matric C in this
**   case has C[1:m, 1:m] diagonal and  C[(m+1):n, 1:n)] is dense. 
**
** arguments are:
**     n         the size of the matrix to be factored
**     m         the size of the diagonal upper portion
**     diag      the diagonal upper portion
**     **matrix  a ragged array containing the dense portion
**     toler     tolerance for detecting singularity
**
**  The diagonal portion of the matrix is unchanged by the factorization.
**  For the dense portion, D occupies the diagonal (of the full matrix).
**  The factorization is returned in the lower triangle.
**   The upper triangle of the matrix is entirely unused by the process (but
**   because of the compressed storage, this isn't much space).
**
**  Return value:  the rank of the matrix (non-negative definite), or -rank
**     if not non-negative definite
**
**  If a column is deemed to be redundant, then that diagonal is set to zero.
**
**   Terry Therneau
*/
#include "survS.h"
#include "survproto.h"

int cholesky3(double **matrix, int n, int m, double *diag, double toler)
    {
    double temp;
    int  i,j,k;
    double eps, pivot;
    int rank;
    int n2;
    int nonneg;

    n2 = n-m;    /* number of full covariates */
   
    nonneg=1;
    eps =0;
    for (i=0; i<m; i++) if (diag[i] <eps) eps = diag[i];
    for (i=0; i<n2; i++) if (matrix[i][i+m] > eps)  eps = matrix[i][i+m];
    eps *= toler;

    rank =0;
    /* pivot out the diagonal elements */
    for (i=0; i<m; i++) {
	pivot = diag[i];
        if (pivot < eps) {
            for (j=0; j<n2; j++) matrix[j][i] =0;
            if (pivot < -8*eps) nonneg= -1;
            }
	else {
	    rank++;
	    for (j=0; j<n2; j++) {
		temp = matrix[j][i] / pivot;
		matrix[j][i] = temp;
		matrix[j][j+m] -= temp*temp*pivot;
		for (k=(j+1); k<n2; k++) matrix[k][j+m] -= temp*matrix[k][i];
		}
	    }
        }

    /* Now the rest of the matrix */
    for (i=0; i<n2; i++) {
	pivot = matrix[i][i+m];
	if (pivot < eps) {
	    for (j=i; j<n2; j++) matrix[j][i+m] =0;  /* zero the column */
            if (pivot < -8*eps) nonneg= -1;
	    }
	else  {
	    rank++;
	    for (j=(i+1); j<n2; j++) {
		temp = matrix[j][i+m]/pivot;
		matrix[j][i+m] = temp;
		matrix[j][j+m] -= temp*temp*pivot;
		for (k=(j+1); k<n2; k++) matrix[k][j+m] -= temp*matrix[k][i+m];
		}
	    }
	}
    return(rank * nonneg);
    }
