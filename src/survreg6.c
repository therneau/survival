/* $Id$ */
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
**	debug   -  >0 causes tracing information.  Can be removed
**	expr    - for callback, the expression to be evaluated
**      rho     - for callback, the environment (R) or frame (Splus) in which
**	            to do the evaluation.  
**  Output
**      beta    - the final coef vector
**      iter    - the number of iterations consumed
**      imat    - the information matrix 
**      loglik  - the final log-liklihood
**      flag    - success flag  0 =ok
**                              -1= did not converge
**      u      - the score vector
**
**  Work arrays
**      newbeta(nvar)- always contains the "next iteration"
**      JJ = the approx variance matrix J'J, guarranteed non-singular
*/
#include "survS.h"
#include "survproto.h"

SEXP survreg6(SEXP maxiter2,   SEXP nvarx,  SEXP y,
	      SEXP ny2,        SEXP covar2, SEXP wtx,
	      SEXP offset2,    SEXP beta2,  SEXP nstratx,
	      SEXP stratax,    SEXP epsx,   SEXP  tolx,
	      SEXP dist,       SEXP dexpr,  SEXP rho) {
    int i,j;	
    int n, maxiter,
	ny;
    double *newbeta;
    int halving, iter;
    double newlk;
    double *loglik, eps, tol_chol;
    double *beta;
    Sint   *flag;
    SEXP   out_beta;
    int    nvar, nvar2, nstrat;
    double **covar;
    Sint   *strat ;
    double *time2, *time1, *status;
    double *offset;
    double **imat, **JJ;
    double *u, *wt, *usave;
    double (*dolik)();   /* will be pointed to survregc1 or survregc2 */
    SEXP  z;
    double *zptr = NULL;
    SEXP  out_iter, out_loglik, out_imat, out_flag;
    SEXP  out_u;
    SEXP  rlist, rlistnames;
    Sint *iter2;
    int nprotect;

    /*
    ** The only input arg that is overwritten is beta
    */
    out_beta = PROTECT(duplicate(beta2));
    beta    = REAL(out_beta);
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

    /*
    ** nvar = # of "real" x variables, for iteration
    ** nvar2= # of parameters to maximize = nvar + nstrat
    ** nstrat= # of strata, where 0== fixed sigma
    */
    nvar2 = nvar + nstrat;   /* number of coefficients */

    /*
    ** Create the output variables
    */
    PROTECT(out_imat = allocVector(REALSXP, nvar2*nvar2)); 
    imat = dmatrix(REAL(out_imat), nvar2, nvar2);
    PROTECT(out_iter   = allocVector(INTSXP, 1));
    iter2 = INTEGER(out_iter);
    PROTECT(out_loglik = allocVector(REALSXP, 1));
    loglik = REAL(out_loglik);
    PROTECT(out_flag   = allocVector(INTSXP, 1));
    flag = INTEGER(out_flag);
    PROTECT(out_u = allocVector(REALSXP, nvar2));
    usave = REAL(out_u);
    nprotect = 6;

    /* Create scratch variables 
    **  u = working version of score vector, overwritten with u H-inv during
    **          Newton steps
    ** usave = a copy of u, after each Newton step.  Returned to the S 
    **   parent routine, and also used to "backtrack" when we need to fail
    **   over to a Fisher step after NR + halving didn't work
    */
    newbeta = (double *) Calloc(LENGTH(beta2) + nvar2 + nvar2*nvar2, double);
    u = newbeta + length(beta2);
    JJ  = dmatrix(u +nvar2, nvar2, nvar2);

    /* 
    ** fixed scale parameters were tacked onto the end of beta at input
    **  copy them to to the end of newbeta as well (survregc1/c2 expects em)
    */
    for (i=nvar; i<LENGTH(beta2); i++) newbeta[i] = beta[i];
       

    if (ny==2) {
	time1= REAL(y); 
	status = time1 +n;
	time2 = NULL;    /*keep gcc from complaining */
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
	**  needs to be of length n + number of status==3 obs
	*/
	j=0;
	for (i=0; i<n; i++) if (status[i]==3) j++;
	PROTECT(z = allocVector(REALSXP, n+j));
	defineVar(install("z"), z, rho);
	zptr = REAL(z);
	nprotect = nprotect+1;
	}
  
    /*
    ** do the initial iteration step
    */
    *loglik = (*dolik)(n,      nvar,             nstrat,  0,
		       beta,   asInteger(dist),  strat,   offset,
		       time1,  time2,            status,  wt,
		       covar,  imat,             JJ,      u,
		       dexpr,  rho,              zptr,
		       0,      NULL,             NULL,    NULL);
    for (i=0; i<nvar2; i++) usave[i] = u[i];
    /*
    ** Why cholesky3 (with 0 sparse terms) instead of cholesky2, which assumes
    **  0 terms?  Because cholesky2 expects the data be in the upper triangle,
    **  cholesky3 expects it in the lower triangle.  (Bad design on my part, 
    **  but long years ago now).  survreg7 also calls survregc1/2, and it does
    **  potentially have sparse terms.  Both choleskys return the L matrix
    **  in the lower triangle, so I can follow with chsolve2 and chinv2.
    */
    *flag= cholesky3(imat, nvar2, 0, NULL, tol_chol);
    if (*flag < 0) {
	i = cholesky3(JJ, nvar2, 0, NULL, tol_chol);
	chsolve2(JJ, nvar2, u);
	}
    else chsolve2(imat,nvar2,u);        /* a replaced by  a *inverse(i) */
   
    /*
    **  Never, never complain about convergence on the first step.  That way,
    **  if someone HAS to they can force one iter at a time.
    */
    for (i=0; i<nvar2; i++) {
	newbeta[i] = beta[i] + u[i];
	}
    if (maxiter==0) {
	chinv2(imat,nvar2);
	for (i=1; i<nvar2; i++)
	    for (j=0; j<i; j++)  imat[i][j] = imat[j][i];
	*iter2 =0;
	goto alldone;
	}


    /*
    ** here is the main loop
    */
    halving =0 ;             /* >0 when in the midst of "step halving" */
    newlk = (*dolik)(n,      nvar,             nstrat,  0,
		     newbeta,asInteger(dist),  strat,   offset,
		     time1,  time2,            status,  wt,
		     covar,  imat,             JJ,      u,
		     dexpr,  rho,              zptr,
		     0,      NULL,             NULL,    NULL); 
    for (i=0; i<nvar2; i++) usave[i] = u[i];

    for (iter=1; iter<= maxiter; iter++) {
	/* 
	**   Am I done?  Check for convergence, then update betas
	*/
	if (fabs(1-(*loglik/newlk))<=eps && halving==0 ) { /* all done */
	    *loglik = newlk;
	    *flag = cholesky3(imat, nvar2, 0, NULL, tol_chol);

	    chinv2(imat, nvar2);     /* invert the information matrix */
	    for (i=1; i<nvar2; i++){
		for (j=0; j<i; j++)  imat[i][j] = imat[j][i];
	        }
	    for (i=0; i<nvar2; i++)
		beta[i] = newbeta[i];
	    if (halving==1) *flag= 1000; /*didn't converge after all */
	    *iter2 = iter;
	    goto alldone;
	    }
	
	if (newlk < *loglik)   {    /*it is not converging ! */
	    for (j=0; j<5 && newlk < *loglik; j++) {
		halving++;
		for (i=0; i<nvar2; i++)
		    newbeta[i] = (newbeta[i] + beta[i]) /2; 
		/*
		** Special code for sigmas.  Often, they are the part
		**  that gets this routine in trouble.  The prior NR step
		**  may have decreased one of them by a factor of >10, in which
		**  case step halving isn't quite enough.  Make sure the new
		**  try differs from the last good one by no more than 1/3
		**  approx log(3) = 1.1
		**  Step halving isn't enough of a "back away" when a
		**  log(sigma) goes from 0.5 to -3, or has become singular.
		*/
		if (halving==1) {  /* only the first time */
		    for (i=0; i<nstrat; i++) {
			if ((beta[nvar+i]-newbeta[nvar+i])> 1.1)
			    newbeta[nvar+i] = beta[nvar+i] - 1.1;  
			}
		    }
		newlk = (*dolik)(n,      nvar,             nstrat,  1,
				 newbeta,asInteger(dist),  strat,   offset,
				 time1,  time2,            status,  wt,
				 covar,  imat,             JJ,      u,
				 dexpr,  rho,              zptr,
				 0,      NULL,             NULL,    NULL); 
		}
	    }

	else {    /* take a standard NR step */
	    halving=0;
	    *loglik = newlk;
	    *flag = cholesky3(imat, nvar2, 0, NULL, tol_chol);
	    if (*flag < 0) {
		i = cholesky3(JJ, nvar2, 0, NULL, tol_chol);
		chsolve2(JJ, nvar2, u);
		}
	    else chsolve2(imat,nvar2,u);
	    for (i=0; i<nvar2; i++) {
		beta[i] = newbeta[i];
		newbeta[i] = newbeta[i] +  u[i];
		}
	    }
	
	newlk = (*dolik)(n,      nvar,             nstrat,  0,
			 newbeta,asInteger(dist),  strat,   offset,
			 time1,  time2,            status,  wt,
			 covar,  imat,             JJ,      u,
			 dexpr,  rho,              zptr,
			 0,      NULL,             NULL,    NULL); 
	for (i=0; i<nvar2; i++) usave[i] = u[i];
	}   /* return for another iteration */
    *iter2 = maxiter;
    *loglik = newlk;
    cholesky3(imat, nvar2, 0, NULL, tol_chol);
    chinv2(imat, nvar2); 
    for (i=1; i<nvar2; i++) {
        for (j=0; j<i; j++)  imat[i][j] = imat[j][i];
        }
    
    for (i=0; i<nvar2; i++)
	beta[i] = newbeta[i];
    *flag= 1000;

    /*
    ** Put together the return list
    */
alldone:
    
    PROTECT(rlist=allocVector(VECSXP, 6));
    SET_VECTOR_ELT(rlist, 0, out_beta);
    SET_VECTOR_ELT(rlist, 1, out_iter);
    SET_VECTOR_ELT(rlist, 2, out_imat);
    SET_VECTOR_ELT(rlist, 3, out_loglik);
    SET_VECTOR_ELT(rlist, 4, out_flag);
    SET_VECTOR_ELT(rlist, 5, out_u);

    /* add names to the objects */
    PROTECT(rlistnames = allocVector(STRSXP, 6));
    SET_STRING_ELT(rlistnames, 0, mkChar("coef"));
    SET_STRING_ELT(rlistnames, 1, mkChar("iter"));
    SET_STRING_ELT(rlistnames, 2, mkChar("var"));
    SET_STRING_ELT(rlistnames, 3, mkChar("loglik"));
    SET_STRING_ELT(rlistnames, 4, mkChar("flag"));
    SET_STRING_ELT(rlistnames, 5, mkChar("u"));
#ifdef USING_R
    setAttrib(rlist, R_NamesSymbol, rlistnames);
#else
    /*
    ** In Splus 8.0.1, Rinternals.h strings don't yet work in all cases 
    */  
    SET_NAMES(rlist, rlistnames);
#endif    

    UNPROTECT(nprotect + 2);
    Free(newbeta);
    return(rlist);
    }
    
