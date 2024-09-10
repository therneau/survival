/*
** General cholesky decompostion
*/
#include "survS.h"
#include "survproto.h"

SEXP gchol(SEXP matrix2, SEXP toler2) {
    int i,j;
    int n;
    double **mat;
    double *toler;
    SEXP gc;   /* the returned matrix */
    
    PROTECT(gc = duplicate(matrix2));
    n = nrows(gc);
    mat = dmatrix(REAL(gc), n, n);
    toler = REAL(toler2);
    
    i = cholesky5(mat, n, *toler);
 
    /* zero out the upper triangle */
    for (i=0; i<n; i++) {
	for (j= i+1; j<n; j++) mat[i][j] =0;
	}
    UNPROTECT(1);
    return(gc);
    }


SEXP gchol_solve(SEXP x2, SEXP y2, SEXP flag2) {
    int n;
    double **mat;
    int flag;
    
    SEXP new_;   /*returned matrix */
    
    n = nrows(x2);
    flag = asInteger(flag2);
    PROTECT(new_ = duplicate(x2));

    mat = dmatrix(REAL(new_), n, n);
    chsolve5(mat, n, REAL(y2), flag);
    UNPROTECT(1);
    return(new_);
    }
    
SEXP gchol_inv(SEXP matrix, SEXP flag2) {
    int n;
    double **mat;
    int i,j;
    int flag;
    SEXP new_;  /* returned matrix */

    n = nrows(matrix);
    flag = asInteger(flag2);
    PROTECT(new_ = duplicate(matrix));
    mat = dmatrix(REAL(new_), n, n);

    chinv5(mat, n, flag);

    /*
    **  the result of chinv5 has the inverse of L and full inverse
    **  all packed together
    */
    if (flag ==1) {
	/* 
	** return L-inverse, by zeroing out the other part
	*/
	for (i=0; i<n; i++) {
	    mat[i][i] = 1;
	    for (j=i+1; j<n; j++) mat[i][j] =0;
	    }
	}
    else {
	/* 
	** replicate the lower part into the upper one, for a symmetric result
	*/
	for (i=0; i<n; i++) {
	    for (j=i+1; j<n; j++) mat[j][i] = mat[i][j];
	    }
	}
    
    UNPROTECT(1);
    return(new_);
    }
   
