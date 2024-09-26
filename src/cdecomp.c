/*
** Compute the eigenvectors for the upper triangular matrix R
*/
#include <math.h>
#include "R.h"
#include "Rinternals.h"

SEXP cdecomp(SEXP R2, SEXP time2) {
    int i,j,k;
    int nc, ii;
    
    static const char *outnames[]= {"d", "A", "Ainv", 
                                    "P", ""};    
    SEXP rval, stemp;
    double *R, *A, *Ainv, *P;
    double *dd, temp, *ediag;
    double time;

    nc = ncols(R2);   /* number of columns */
    R = REAL(R2);
    time = asReal(time2);

    /* Make the output matrices as copies of R, so as to inherit
    **   the dimnames and etc
    */
    
    PROTECT(rval = mkNamed(VECSXP, outnames));
    stemp=  SET_VECTOR_ELT(rval, 0, allocVector(REALSXP, nc));
    dd = REAL(stemp);
    stemp = SET_VECTOR_ELT(rval, 1, allocMatrix(REALSXP, nc, nc));
    A = REAL(stemp);
    for (i =0; i< nc*nc; i++) A[i] =0;   /* R does not zero memory */
    stemp = SET_VECTOR_ELT(rval, 2, duplicate(stemp));
    Ainv = REAL(stemp);
    stemp = SET_VECTOR_ELT(rval, 3, duplicate(stemp));
    P = REAL(stemp);
   
    ediag = (double *) R_alloc(nc, sizeof(double));
    
    /* 
    **        Compute the eigenvectors
    **   For each column of R, find x such that Rx = kx
    **   The eigenvalue k is R[i,i], x is a column of A
    **  Remember that R is in column order, so the i,j element is in
    **   location i + j*nc
    */
    ii =0; /* contains i * nc */
    for (i=0; i<nc; i++) { /* computations for column i */
        dd[i] = R[i +ii];    /* the i,i diagonal element = eigenvalue*/
        A[i +ii] = 1.0;
        for (j=(i-1); j >=0; j--) {  /* fill in the rest */
            temp =0;
            for (k=j; k<=i; k++) temp += R[j + k*nc]* A[k +ii];
            A[j +ii] = temp/(dd[i]- R[j + j*nc]);
        }
        ii += nc;
    }
    
    /*
    ** Solve for A-inverse, which is also upper triangular. The diagonal
    **  of A and the diagonal of A-inverse are both 1.  At the same time 
    **  solve for P = A D Ainverse, where D is a diagonal matrix 
    **  with exp(eigenvalues) on the diagonal.
    ** P will also be upper triangular, and we can solve for it using
    **  nearly the same code as above.  The prior block had RA = x with A the
    **  unknown and x successive colums of the identity matrix. 
    **  We have PA = AD, so x is successively columns of AD.
    ** Imagine P and A are 4x4 and we are solving for the second row
    **  of P.  Remember that P[2,1]= A[2,3] = A[2,4] =0; the equations for
    **  this row of P are:
    **
    **    0*A[1,2] + P[2,2]A[2,2] + P[2,3] 0     + P[2,4] 0     = A[2,2] D[2]
    **    0*A[1,3] + P[2,2]A[2,3] + P[2,3]A[3,3] + P[2,4] 0     = A[2,3] D[3]
    **    0*A[1,4] + P[2,2]A[2,4] + P[2,3]A[3,4] + P[2,4]A[4,4] = A[2,4] D[4]
    **
    **  For A-inverse the equations are (use U= A-inverse for a moment)
    **    0*A[1,2] + U[2,2]A[2,2] + U[2,3] 0     + U[2,4] 0     = 1
    **    0*A[1,3] + U[2,2]A[2,3] + U[2,3]A[3,3] + U[2,4] 0     = 0
    **    0*A[1,4] + U[2,2]A[2,4] + U[2,3]A[3,4] + U[2,4]A[4,4] = 0
    */
    
    ii =0; /* contains i * nc */
    for (i=0; i<nc; i++) ediag[i] = exp(time* dd[i]);
    for (i=0; i<nc; i++) { 
        /* computations for column i of A-inverse */
        Ainv[i+ii] = 1.0 ;
        for (j=(i-1); j >=0; j--) {  /* fill in the rest of the column*/
            temp =0;
            for (k=j+1; k<=i; k++) temp += A[j + k*nc]* Ainv[k +ii];
            Ainv[j +ii] = -temp;
        }
        
        /* column i of P */
        P[i + ii] = ediag[i];
        for (j=0; j<i; j++) {
            temp =0;
            for (k=j; k<nc; k++) temp += A[j + k*nc] * Ainv[k+ii] * ediag[k];
            P[j+ii] = temp;
        }
        
        /* alternate computations for row i of P, does not use Ainv*/
        /*P[i +ii] = ediag[i];
          for (j=i+1; j<nc; j++) { 
              temp =0;
              for (k=i; k<j; k++) temp += P[i+ k*nc]* A[k + j*nc];
              P[i + j*nc] = (A[i + j*nc]*ediag[j] - temp)/A[j + j*nc];
          } 
        */
        ii += nc;
    }
    UNPROTECT(1);
    return(rval);
}
