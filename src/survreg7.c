/*
** Fit one of several censored data distributions
**
** Input
**      maxiter - max # of iterations allowed
**      n       - number of subjects
**      nvar    - number of variables in the x matrix
**      y       - matrix of start, stop, event
**      ny      - # columns of y.  =3 if there is interval censored data
**      event   - 1=exact, 0= right censored, 2=left censored, 3=interval
**      covar   - covariates for patient i
**                   Note that S sends this in column major order
**      wt      - vector of case weights (usually 1)
**      offset  - offset vector (usually 0)
**      beta    - initial values for the parameters, of length 1+nvar
**		     last element contains the scale
**      nstrat  - an indicator: 0= scale is fixed, 1= estimate scale,
**		    >1: estimate multiple scales (strata)
**      strat   - if nstrat>0, contains the strata number for each subject
**      eps     - tolerance for convergence.  Iteration continues until the
**                  relative change in the deviance is <= eps.
**      tol_chol- tolerance for Cholesky decomposition
**      dist    -  1=extreme value, 2=logistic, 3=gaussian, 4=callback
**	dexpr   - for callback, the expression to be evaluated that evaluates
**                  the distribution function of the random effect
**      rho     - for callback, the environment (R) or frame (Splus) in which
**	            to do evaluations 
**      ptype   - 1= sparse penalties, 2=dense penalties, 1+2 = both, 0=none
**      pdiag   - 0 = the penalty matrix is diagonal
**      nfrail  - number of levels of the sparse term (0 = no sparse term)
**      fgrp    - which frailty group each subject is in
**      pexpr1  - for callback, the expression to eval for sparse penalties
**      pexpr2  - the expression for dense penalties
**      
**  Output
**      beta    - the final coef vector
**      iter    - the number of iterations consumed
**      hmat    - the cholesky of the penalized information matrix 
**      hinv    - the cholesky of the inverse of hmat
**      hdiag   - diagonal portion of hinv
**      loglik  - the final log-liklihood
**      u       - the final score vector.  Usually =0 at convergence, but
**                  useful in other cases for a score test.
**      flag    - success flag  0 =ok
**                              -1= did not converge
**
**  Work arrays
**      newbeta(nvar)- always contains the "next iteration"
**      u(nvar)      - first deriv of the loglik
**      JJ = the approx variance matrix J'J, guarranteed non-singular
**
**  Notes on hmat:  H will be p=(nfrail+ nvar + nstrat) square, but the
**    upper left nfrail*nfrail corner is a diagonal matrix.  It is stored
**    as "hdiag", which "hmat" contains the remaining dense portion. If 
**    H = LDL' (see cholesky3), then H-inverse = (L-inv)' (Dinv) (L-inv)
**    where D is diagonal and L is lower-triangular with ones on the diagonal,
**    and L[1:nfrail, 1:nfrail] is the identity.
**    The return parts are hmat = L[(nfrail+1):p, 1:p], hinv L-inverse[ same],
**    and g=hdiag = D-inverse.  See coxpenal.df for more.
*/
#include "survS.h"
#include "survproto.h"

SEXP survreg7(SEXP maxiter2,   SEXP nvarx,  SEXP y,
	      SEXP ny2,        SEXP covar2, SEXP wtx,
	      SEXP offset2,    SEXP beta2,  SEXP nstratx,
	      SEXP stratax,    SEXP epsx,   SEXP tolx,
	      SEXP dist,       SEXP dexpr,  SEXP rho,
	      SEXP ptype2,     SEXP pdiag2, SEXP nfrail2,
	      SEXP fgrp2,      SEXP pexpr1, SEXP pexpr2) {
    /* local variables */
    int i,j;	
    int nvar, nvar2, nvar3, nstrat;
    int iter;
    double newlk =0;
    double (*dolik)();   /* will point to (*dolik) or survregc2 */
    double x1, x2, x3, x4;
    double y1, y2, y3;
    int golden, goright;
    double newpen;

    /* pointers for the data regions of the input arguments */
    double **covar;
    Sint   *strat ;
    double *time2, *time1, *status;
    double *offset;
    Sint *fgrp;
    double *wt;

    /* copies of the scalar input arguments */
    double eps, tol_chol;
    int n, maxiter, ny;
    int nfrail, ptype, pdiag;

    /* Variables allocated in this routine */
    double *jdiag, 
	   *newbeta,
	   *u;

    /* variables for the callback code */
    SEXP   coef1,  coef2;
    double *cptr1=NULL, *cptr2 =NULL;    /* stop a gcc warning */
    double **JJ;
    SEXP  z;
    double *zptr = NULL;

    /* structures and pointers for the returned list object */
    SEXP  out_iter, out_loglik, out_hmat, 
	  out_hinv, out_flag,   out_beta;
    SEXP  out_penalty;
    SEXP  out_hdiag, out_u;
    double *loglik, *usave;
    double **hmat, **hinv, *beta, *hdiag;
    double *penalty;
    SEXP  rlist, rlistnames;
    Sint *iter2, *flag;
    int nprotect;   /* number of PROTECT calls that I have issued */

    /*
    ** The only input arg that is rewritten is beta, so no need to duplicate
    */
    maxiter = asInteger(maxiter2);
    n =       LENGTH(wtx);
    ny =      asInteger(ny2);
    nvar =    asInteger(nvarx);
    offset =  REAL(offset2);
    nstrat =  asInteger(nstratx);
    strat  =  INTEGER(stratax);
    wt =      REAL(wtx);
    eps    =  asReal(epsx);
    tol_chol= asReal(tolx);
    covar =   dmatrix(REAL(covar2), n, nvar);
    nfrail =  asInteger(nfrail2);
    ptype  =  asInteger(ptype2);
    pdiag  =  asInteger(pdiag2);
    fgrp   =  INTEGER(fgrp2);

    /*
    ** nvar = # of "real" x variables, found in the coefficient matrix
    ** nvar2= size of the dense portion of hmat = nvar + nstrat
    ** nvar3= #coefficients = nfrail + nvar2
    ** nstrat= # of strata, where 0== fixed sigma
    */
    nvar2 = nvar + nstrat;   /* number of coefficients */
    nvar3 = nvar2 + nfrail;

    /*
    ** Create the output variables
    */
    PROTECT(out_beta = duplicate(beta2));
    beta    = REAL(out_beta);
    PROTECT(out_hmat = allocVector(REALSXP, nvar3*nvar2)); 
    hmat = dmatrix(REAL(out_hmat), nvar3, nvar2);
    PROTECT(out_hinv = allocVector(REALSXP, nvar3*nvar2)); 
    hinv = dmatrix(REAL(out_hinv), nvar3, nvar2);
    PROTECT(out_hdiag = allocVector(REALSXP, nvar3));
    hdiag = REAL(out_hdiag);
    PROTECT(out_iter   = allocVector(INTSXP, 1));
    iter2 = INTEGER(out_iter);
    PROTECT(out_loglik = allocVector(REALSXP, 1));
    loglik = REAL(out_loglik);
    PROTECT(out_flag   = allocVector(INTSXP, 1));
    flag = INTEGER(out_flag);
    PROTECT(out_u = allocVector(REALSXP, nvar3));
    usave = REAL(out_u);  /* the working vector 'u' gets destroyed in chsolve*/
    PROTECT(out_penalty= allocVector(REALSXP, 1));
    penalty = REAL(out_penalty);
    nprotect =9;

    /* Create the scratch vectors 
    **  u = working version of score vector, overwritten with u H-inv during
    **          Newton steps
    ** usave = a copy of u, after each Newton step.  Returned to the S 
    **   parent routine, and also used to "backtrack" when we need to fail
    **   over to a Fisher step instead of an NR step
    */
    newbeta = Calloc(LENGTH(beta2) + nvar3 + nfrail + nvar2*nvar3, double);
    jdiag = newbeta + length(beta2);
    u  = jdiag + nfrail;
    JJ  = dmatrix(u + nvar3,  nvar3, nvar2);

    /* 
    ** fixed scale parameters were tacked onto the end of beta at input
    **  copy them to to the end of newbeta as well ((*dolik) expects them)
    */
    for (i=nvar; i<LENGTH(beta2); i++) newbeta[i] = beta[i];

    if (ny==2) {
	time1= REAL(y);
	status = time1 +n;
        time2 = NULL;       /*quiet a compiler warning*/
	}
    else {
	time1= REAL(y);
	time2 = time1 + n;
	status = time2 +n;
	}

    if (asInteger(dist) <4) {
	/* case 1, the "built in" distributions */
	dolik = survregc1;
	}
    else {
        dolik = survregc2;
	/*
	** Create the vector z in the contained data frame
	**  (dexpr is implicitly a function of "z")
	** Technically it is not "in" rho, just "visible from" rho
	*/
	PROTECT(z = allocVector(REALSXP, n));
	defineVar(install("z"), z, rho);
	zptr = REAL(z);
	nprotect++;
	}

    if (ptype==1 || ptype==3) {
	/* 
	**  There is a penalty on the sparse terms
	**  Create the vector coef1 in the contained data frame
	**  (pexpr1 is implicitly a function of 'coef1')
	*/
	PROTECT(coef1 = allocVector(REALSXP, nfrail));
	defineVar(install("coef1"), coef1, rho);
	cptr1 = REAL(coef1);
	nprotect++;
	}
    if (ptype >1) {
	/*
	** There is a penalty on the non-sparse terms
	** Create the vector coef2 in the contained data frame
	** (pexpr2 is implicitly a function of 'coef2')
	** Since scale parameters are never panalized, only the first
	**   nvar of the coefficients are passed in.
	*/
	PROTECT(coef2 = allocVector(REALSXP, nvar));
	defineVar(install("coef2"), coef2, rho);
	cptr2 = REAL(coef2);
	nprotect++;
	}

    /*
    ** Get the loglik, score, and hessian for the initial parameters
    */
    *loglik = (*dolik)(n,      nvar,             nstrat,  0,
		       beta,   asInteger(dist),  strat,   offset,
		       time1,  time2,            status,  wt,
		       covar,  hmat,             JJ,      u,
		       dexpr,  rho,              zptr,
		       nfrail, fgrp,             hdiag,   jdiag);
    survpenal(0, nfrail, nvar,  hmat,  JJ,     hdiag, jdiag,  u,     beta,  
	         penalty, ptype, pdiag, pexpr1, cptr1, pexpr2, cptr2, rho);
    *loglik += *penalty;
    for (i=0; i<nvar3; i++) usave[i] = u[i];

    /*
    ** The code below is careful to return a consistent answer: beta, 
    **   u(beta) and hmat(beta); e.g., not newbeta and u(beta)
    **   This isn't an issue when the algorithm converges, since by
    **   definition old beta = new beta then.  But it is an issue for
    **   certain cases where a user wants exactly 0 or 1 iteration.
    */  
    for (iter=1; iter <=maxiter; iter++) {
	/*
	** Take a NR or Fisher step to get newbeta
	*/
	*flag = cholesky3(hmat, nvar3, nfrail, hdiag, tol_chol);
	if (*flag <0) {
	    /* Fisher step */
	    cholesky3(JJ, nvar3, nfrail, jdiag, tol_chol);
	    chsolve3(JJ, nvar3, nfrail, jdiag, u);
	    }
	else {  /* Newton-Raphson step */
	    chsolve3(hmat,nvar3, nfrail, hdiag, u);
	    }

	for (i=0; i<nvar3; i++) {
	    newbeta[i] = beta[i] + u[i];
	    }

	newlk = (*dolik)(n,       nvar,             nstrat,  0,
			  newbeta, asInteger(dist),  strat,   offset,
			  time1,   time2,            status,  wt,
			  covar,   hmat,             JJ,      u,
			  dexpr,   rho,              zptr,
			  nfrail,  fgrp,             hdiag,   jdiag);
	survpenal(0, nfrail, nvar,  hmat,  JJ,     hdiag, jdiag,  u, newbeta,  
		  &newpen, ptype, pdiag, pexpr1, cptr1, pexpr2, cptr2, rho);
	newlk += newpen;

	/* 
	**   Am I done?  If so we can drop out of the iter loop
	*/
	if (fabs(1-(*loglik/newlk))<=eps) { /* all done */
	    *loglik = newlk;
	    *penalty= newpen;
	    for (i=0; i<nvar3; i++) {
		beta[i] = newbeta[i];
		usave[i]= u[i];
		}
	    *iter2 = iter;
	    goto alldone;
	    }
	
	/*
	**  Make sure that the step was an improvement
	*/
	if (newlk < *loglik)   {    /*it is not converging ! */
	    /*
	    ** Do 8 steps of a golden section search - just a little more
	    **    sophisticated than step halving.  We optimize over
	    **    newbeta = coef + alpha * u, for the best alpha.  Since
	    **    this is just an intermediate correction, there is no
	    **    need for the exact optimum.
	    ** We have found that the sigmas are often the
	    **   cause of trouble: the prior NR step may have decreased one
	    **   of them by a factor of up to 10, leading to a yhat that is
	    **   a huge outlier in the density.  So we first modify u to
	    **   ensure that no sigma shrinkage is more than about 1/3, i.e.,
	    **   a -1.1 change in log sigma.  
	    ** We use "u" as a scratch var, holding the last good step.  Note
	    **   the use of whichcase=1 until we're done, so u is left alone.
	    ** We sometimes see saddle points where alpha <0, btw.
	    */
	    for (i=0; i<nvar3; i++) u[i] = newbeta[i] - beta[i];
	    for (i=0; i<nstrat; i++) 
		if (u[i+nvar+nfrail] < -0.7) u[i+nvar+nfrail] = -1.1;

	    /*
  	    ** Setup: x1=left endpoint of search, x2,x3= middle, x4=right
	    **   First we bracket the root
	    */
	    x4=1;
	    x1 = 0;
	    x3 = 1;
	    y1 = *loglik;
	    y3 = newlk;
	    while (y1 > y3) {
		x4 = x3;
		x3 = x1;
		x1 = x3 - (x4-x3)/.618;
		y3 = y1;
		for (i=0; i<nvar3; i++)
		    newbeta[i] = beta[i] + u[i]*x1;

		y1 = (*dolik)(n,      nvar,             nstrat,  1,
			      newbeta,asInteger(dist),  strat,   offset,
			      time1,  time2,            status,  wt,
			      covar,  hmat,             JJ,      u,
			      dexpr,  rho,              zptr,
			      nfrail, fgrp,             hdiag,   jdiag);
		survpenal(1, nfrail,  nvar,   hmat,  JJ, hdiag, jdiag,  
			  u, newbeta, &newpen, ptype, pdiag,
			  pexpr1, cptr1, pexpr2, cptr2, rho);
		y1 += newpen;
		}
	    
	    x2 = .618*x3 + .382*x1;
	    for (i=0; i<nvar3; i++)
		    newbeta[i] = beta[i] + u[i]*x2;
	    y2 = (*dolik)(n,      nvar,             nstrat,  1,
			  newbeta,asInteger(dist),  strat,   offset,
			  time1,  time2,            status,  wt,
			  covar,  hmat,             JJ,      u,
			  dexpr,  rho,              zptr,
			  nfrail, fgrp,             hdiag,   jdiag);
	    survpenal(1, nfrail,  nvar,   hmat,  JJ, hdiag, jdiag,  
		      u, newbeta, &newpen, ptype, pdiag,
		      pexpr1, cptr1, pexpr2, cptr2, rho);
	    y2 += newpen;

	    for (golden=0; golden< 8; golden++) {
		if (y3 > y2) { /* toss away the interval from x1 to x2 */
		    x1=x2;
		    x2=x3; 
		    x3 = .618*x4 + .382*x1;
		    y2 =y3;
		    for (i=0; i<nvar3; i++)
			newbeta[i] = beta[i] + u[i]*x3;
		    goright=1;
		    }
		else { /* toss away the interval from x3 to x4 */
		    x4 = x3;
		    x3 = x2;
		    x2 = .618*x1 + .382*x4;
		    y3 = y2;
		    for (i=0; i<nvar3; i++)
			newbeta[i] = beta[i] + u[i]*x2;
		    goright =0;
		    }
		     
		newlk = (*dolik)(n,      nvar,             nstrat,  1,
				 newbeta,asInteger(dist),  strat,   offset,
				 time1,  time2,            status,  wt,
				 covar,  hmat,             JJ,      u,
				 dexpr,  rho,              zptr,
				 nfrail, fgrp,             hdiag,   jdiag);
		survpenal(1, nfrail,  nvar,   hmat,  JJ, hdiag, jdiag,  
			  u, newbeta, &newpen, ptype, pdiag,
			  pexpr1, cptr1, pexpr2, cptr2, rho);
		     
		if (goright) y3= newlk + newpen;
		else         y2= newlk + newpen;
	        }

	     if (y2 > *loglik || y3 > *loglik) {
		 /* Success - keep the better guess & compute derivatives */
		 if (y2 > y3) {
		     for (i=0; i<nvar3; i++) newbeta[i] = beta[i] + u[i]*x2;
		     }
		 else {
		     for (i=0; i<nvar3; i++) newbeta[i] = beta[i] + u[i]*x3;
		     }
		 newlk = (*dolik)(n,      nvar,             nstrat,  0,
				  newbeta,asInteger(dist),  strat,   offset,
				  time1,  time2,            status,  wt,
				  covar,  hmat,             JJ,      u,
				  dexpr,  rho,              zptr,
				  nfrail, fgrp,             hdiag,   jdiag);
		 survpenal(0, nfrail,  nvar,   hmat,  JJ, hdiag, jdiag,  
			   u, newbeta, &newpen, ptype, pdiag,
			   pexpr1, cptr1, pexpr2, cptr2, rho);
		 newlk += newpen;

		 }
	     else { /* abject failure */
		 break;  /* fall out of the iteration loop */
		 }
	     }
	
	    
	/*
	** We have a "newbeta" that is an improvement
	**  Keep it.
	*/
	for (i=0; i<nvar3; i++) {
		beta[i] = newbeta[i];
		usave[i] = u[i];
		}
	*loglik = newlk;
	*penalty= newpen;
	}   /* return for another iteration */

    if (maxiter > 1) *flag= 1000;  /* no "non convergence" for 0 or 1 iter */
    *iter2 = iter;

    /*
    ** Put together the return list
    */
alldone:
    *flag = cholesky3(hmat, nvar3, nfrail, hdiag, tol_chol);
    for (i=0; i<nvar2; i++) {
	for (j=0; j<nvar3; j++)  hinv[i][j] = hmat[i][j];
	}
    chinv3(hinv, nvar3, nfrail, hdiag);

    for (i=nfrail; i<nvar3; i++) {       /*nicer output for S user */
	hdiag[i] = hinv[i-nfrail][i];
	hmat[i-nfrail][i] =1;
	hinv[i-nfrail][i] =1;
	for (j=i+1; j<nvar3; j++) {
	    hmat[i-nfrail][j] = 0;
	    hinv[i-nfrail][j] = 0;
	    }
	}
    Free(newbeta);

    /* Create the list object for return */
    PROTECT(rlist=allocVector(VECSXP, 9));
    SET_VECTOR_ELT(rlist, 0, out_beta);
    SET_VECTOR_ELT(rlist, 1, out_iter);
    SET_VECTOR_ELT(rlist, 2, out_hmat);
    SET_VECTOR_ELT(rlist, 3, out_hinv);
    SET_VECTOR_ELT(rlist, 4, out_loglik);
    SET_VECTOR_ELT(rlist, 5, out_flag);
    SET_VECTOR_ELT(rlist, 6, out_hdiag);
    SET_VECTOR_ELT(rlist, 7, out_u);
    SET_VECTOR_ELT(rlist, 8, out_penalty);

    /* add names to the objects */
    PROTECT(rlistnames = allocVector(STRSXP, 9));
    SET_STRING_ELT(rlistnames, 0, mkChar("coef"));
    SET_STRING_ELT(rlistnames, 1, mkChar("iter"));
    SET_STRING_ELT(rlistnames, 2, mkChar("hmat"));
    SET_STRING_ELT(rlistnames, 3, mkChar("hinv"));
    SET_STRING_ELT(rlistnames, 4, mkChar("loglik"));
    SET_STRING_ELT(rlistnames, 5, mkChar("flag"));
    SET_STRING_ELT(rlistnames, 6, mkChar("hdiag"));
    SET_STRING_ELT(rlistnames, 7, mkChar("u"));
    SET_STRING_ELT(rlistnames, 8, mkChar("penalty"));
#ifdef USING_R
    setAttrib(rlist, R_NamesSymbol, rlistnames);
#else
    /*
    ** In Splus 8.0.1, Rinternals.h strings don't yet work in this case,
    **   so use the Splus style call.
    */  
    SET_NAMES(rlist, rlistnames);
#endif    

    UNPROTECT(nprotect + 2);
    return(rlist);
    }
    
