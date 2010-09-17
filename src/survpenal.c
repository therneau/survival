/*
**  Compute the penaltys to the loglik, first, and second derivatives.
**  Then add them on to the relevant terms.
**  It is possible to have one penalty on the sparse terms and another on
**    the dense portion.
**  If ptype =0 this routine is essentially a no-op.
**
**  When whichcase=1, we only update the loglik.  This happens during 
**   step halving in the parent routine, and we want to leave the other
**   matrices alone in case that routine needs to go back and reuse an
**   old result.
*/
#include "survS.h"

void survpenal(int whichcase, int nfrail,    int  nvar,    double **hmat, 
	       double **JJ,   double *hdiag, double *jdiag,
	       double *u,     double *beta,  double *penalty,
	       int ptype,     int pdiag,     SEXP pexpr1,   double *cptr1, 
	       SEXP pexpr2,   double *cptr2, SEXP rho) {

    int i,j, k;
    Sint *flag;
    SEXP plist;
    double *dptr;

    *penalty =0;
    if (ptype ==1 || ptype==3) {
	/*
	** Get the penalty for the sparse term
	*/
	for (i=0; i<nfrail; i++) cptr1[i] = beta[i];  /* set "coef1" in rho */
	PROTECT(plist = eval(pexpr1, rho));
	/*
	** plist is a list, with elements of a new coef, 
	** penalty to the first derivative, penalty to the second,
	** penalty to the loglik, and a logical flag, in that order.
	** Add these into the relevant objects: u, hmat & JJ, loglik.
	*/
	*penalty += asReal(VECTOR_ELT(plist,3));
	if (whichcase==0) {
	    dptr = REAL(VECTOR_ELT(plist,0));
	    for (i=0; i<nfrail; i++) beta[i] = dptr[i];

	    flag =   LOGICAL(VECTOR_ELT(plist, 4));

	    if (flag[0] >0) {
		/*
		** This happens when a penalty is infinite.  (Which itself
		**  is often how a penalty routine signals that it was given
		**  an illegal value for a parameter).  In this case the
		**  updated value of beta will have already been set to 0
		**  via cptr1 above, which is of course the correct solution,
		**  but the U and H use for the parent routine's Newton-Raphson
		**  step are infinite as well.  We force the u and
		**  hmat matrices to dummy values that will cause no update
		**  (none needed) and more importantly no infinte/infinite
		**  arithmetic errors: u=0 and H = identity.  (Only the
		**  relevant columns of each, of course).
		*/
		for (i=0; i<nfrail; i++) {
		    hdiag[i] =1;
		    jdiag[i] =1;
		    u[i] =0;
		    for (j=0; j<nvar; j++) hmat[j][i] =0; 
		    }
	        }
	    else {
		dptr = REAL(VECTOR_ELT(plist, 1));
		for (i=0; i<nfrail; i++) u[i] += dptr[i];
	    
		dptr = REAL(VECTOR_ELT(plist,2));
		for (i=0; i<nfrail; i++) {
		    hdiag[i] += dptr[i];
		    jdiag[i] += dptr[i];
		    }
	        }
	    }
	UNPROTECT(1);
	}

    if (ptype > 1) {
	/*
	** Get the penalty for the dense part of the matrix
	** Note that penalties never apply to the variance terms,
	**   which means that indices go to nvar, not nvar2
	*/
	for (i=0; i<nvar; i++) cptr2[i] = beta[i+nfrail]; /* set fcn var */
	PROTECT(plist = eval(pexpr2, rho));     /* "call" the function */
	*penalty += asReal(VECTOR_ELT(plist,3));  /* loglik penalty */

	if (whichcase==0) {
	    dptr = REAL(VECTOR_ELT(plist,0));   /* updated coefficients */
	    for (i=0; i<nvar; i++) beta[i+nfrail] = dptr[i];

	    dptr = REAL(VECTOR_ELT(plist, 1));       /* first deriv penalty */
	    for (i=0; i<nvar; i++) u[i+nfrail] += dptr[i];
	    
	    dptr = REAL(VECTOR_ELT(plist,2));
	    if (pdiag==0) {
		/* diagonal penalty */
		for (i=0; i<nvar; i++) {
		    JJ[i][i+nfrail] += dptr[i];
		    hmat[i][i+nfrail] += dptr[i];
		    }
	        }
	    else {
		/* full penalty matrix */
		k=0;
		for (i=0; i<nvar; i++) {
		    for (j=nfrail; j<nvar+nfrail; j++) {
			JJ[i][j] += dptr[k];
			hmat[i][j] += dptr[k];
			k++;
		        }
		    }
	        }

	    flag = LOGICAL(VECTOR_ELT(plist, 4));
	    for (i=0; i<nvar; i++) {
		if (flag[i] ==1) {  /* See comments above on flag[0]*/
		    u[i+nfrail] =0;
		    hmat[i][i+nfrail] =1;
		    for (j=0; j<i; j++) hmat[i][j+nfrail] =0;
		    }
	        }
	    }
	UNPROTECT(1);
	}
	
    }
		
