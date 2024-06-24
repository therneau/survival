/*
** Cox regression fit, replacement for coxfit2 in order
**   to be more frugal about memory: specificly that we 
**   don't make copies of the input data.
** 6/2019 : change variable name "time" to "xtime", Sun OS reserves 'time'
**
**  the input parameters are
**
**       maxiter      :number of iterations
**       time(n)      :time of event or censoring for person i
**       status(n)    :status for the ith person    1=dead , 0=censored
**       covar(nv,n)  :covariates for person i.
**                        Note that S sends this in column major order.
**       strata(n)    :marks the strata.  Will be 1 if this person is the
**                       last one in a strata.  If there are no strata, the
**                       vector can be identically zero, since the nth person's
**                       value is always assumed to be = to 1.
**       offset(n)    :offset for the linear predictor
**       weights(n)   :case weights
**       init         :initial estimate for the coefficients
**       eps          :tolerance for convergence.  Iteration continues until
**                       the percent change in loglikelihood is <= eps.
**       chol_tol     : tolerance for the Cholesky decompostion
**       method       : 0=Breslow, 1=Efron
**       doscale      : center and scale each col of X
**
**  returned parameters
**       means(nv)    : vector of column means of X
**       beta(nv)     :the vector of answers (at start contains initial est)
**       u(nv)        :score vector
**       imat(nv,nv)  :the variance matrix at beta=final
**                      (returned as a vector)
**       loglik(2)    :loglik at beta=initial values, at beta=final
**       sctest       :the score test at beta=initial
*       flag         :success flag  1000  did not converge
**                                   1 to nvar: rank of the solution
**       iter         :actual number of iterations used
**
**  work arrays
**       mark(n)
**       wtave(n)
**       a(nvar), a2(nvar)
**       cmat(nvar,nvar)       ragged array
**       cmat2(nvar,nvar)
**       newbeta(nvar)         always contains the "next iteration"
**
**  calls functions:  cholesky2, chsolve2, chinv2
**
**  the data must be sorted by ascending time within strata
*/
#include <math.h>
#include "survS.h"
#include "survproto.h"

/* 
** these arrays are shared with the subroutine at the bottom,
**  but remain unknown outside this source file
*/

static double *a, *a2, **cmat2;
static double *xtime, *weights, *offset;
static int *status, *strata;
static double *u; 
static double **covar, **cmat, **imat;  /*ragged arrays */
static double **null_score; /*score for each patient and each var*/
static double **null_imat; /*hessian matrix under the null*/

static double coxfit6_iter(int nvar, int nused, int method, double *beta, int update_score);

SEXP coxfit6(SEXP maxiter2,  SEXP time2,   SEXP status2, 
	     SEXP covar2,    SEXP offset2, SEXP weights2,
	     SEXP strata2,   SEXP method2, SEXP eps2, 
	     SEXP toler2,    SEXP ibeta,    SEXP doscale2) {

    int i,j, person;
    double temp, temp2;
    double *newbeta, *scale;
    double halving =0, newlk;
    int notfinite;

    /* copies of scalar input arguments */
    int     nused, nvar, maxiter;
    int     method;
    double  eps, toler;
    int *doscale;
   
    /* returned objects */
    SEXP imat2, means2, beta2, u2, loglik2;
    SEXP null_score2, null_imat2;
    double *beta, *loglik, *means;
    SEXP sctest2, flag2, iter2;
    double *sctest;
    int *flag, *iter;
    SEXP rlist, rlistnames;
    int nprotect;  /* number of protect calls I have issued */

    /* get local copies of some input args */
    nused = LENGTH(offset2);
    nvar  = ncols(covar2);
    method = asInteger(method2);
    maxiter = asInteger(maxiter2);
    eps  = asReal(eps2);     /* convergence criteria */
    toler = asReal(toler2);  /* tolerance for cholesky */
    doscale = INTEGER(doscale2);

    xtime = REAL(time2);
    weights = REAL(weights2);
    offset= REAL(offset2);
    status = INTEGER(status2);
    strata = INTEGER(strata2);
    
    /*
    **  Set up the ragged arrays and scratch space
    **  Normally covar2 does not need to be duplicated, even though
    **  we are going to modify it, due to the way this routine was
    **  was called.  But check.
    */
    nprotect =0;
    if (MAYBE_REFERENCED(covar2)) {
	PROTECT(covar2 = duplicate(covar2)); 
	nprotect++;
	}
    covar= dmatrix(REAL(covar2), nused, nvar);

    PROTECT(imat2 = allocVector(REALSXP, nvar*nvar)); 
    nprotect++;
    imat = dmatrix(REAL(imat2),  nvar, nvar);
    a = (double *) R_alloc(2*nvar*nvar + 4*nvar, sizeof(double));
    newbeta = a + nvar;
    a2 = newbeta + nvar;
    scale = a2 + nvar;
    cmat = dmatrix(scale + nvar,   nvar, nvar);
    cmat2= dmatrix(scale + nvar +nvar*nvar, nvar, nvar);

    /* 
    ** create output variables
    */ 
    PROTECT(beta2 = duplicate(ibeta));
    beta = REAL(beta2);
    PROTECT(means2 = allocVector(REALSXP, nvar));
    means = REAL(means2);
    PROTECT(u2 = allocVector(REALSXP, nvar));
    u = REAL(u2);
    PROTECT(loglik2 = allocVector(REALSXP, 2)); 
    loglik = REAL(loglik2);
    PROTECT(sctest2 = allocVector(REALSXP, 1));
    sctest = REAL(sctest2);
    PROTECT(flag2 = allocVector(INTSXP, 1));
    flag = INTEGER(flag2);
    PROTECT(iter2 = allocVector(INTSXP, 1));
    iter = INTEGER(iter2);
    nprotect += 7;
    
    PROTECT(null_score2 = allocVector(REALSXP, nused*nvar)); 
    null_score = dmatrix(REAL(null_score2), nused, nvar);
    nprotect++;
    
    PROTECT(null_imat2 = allocVector(REALSXP, nvar*nvar)); 
    nprotect++;
    null_imat = dmatrix(REAL(null_imat2),  nvar, nvar);
    
    /*
    ** Subtract the mean from each covar, as this makes the regression
    **  much more stable.
    */
    temp2=0;
    for (i=0; i<nused; i++) {
	temp2 += weights[i];
    }	
    for (i=0; i<nvar; i++) {
	if (doscale[i]==0) {scale[i] = 1.0; means[i] =0;}
	else {
	    temp=0;
	    for (person=0; person<nused; person++) 
		temp += weights[person] * covar[i][person];
	    temp /= temp2;
	    means[i] = temp;
	    for (person=0; person<nused; person++) covar[i][person] -=temp;

	    temp =0;
	    for (person=0; person<nused; person++) {
		temp += weights[person] * fabs(covar[i][person]);
	    }
	    if (temp > 0) temp = temp2/temp;   /* scaling */
	    else temp=1.0; /* rare case of a constant covariate */
	    scale[i] = temp;
	    for (person=0; person<nused; person++) {
		covar[i][person] *= temp;
		}
	    }
	}
 
    for (i=0; i<nvar; i++) beta[i] /= scale[i]; /*rescale initial betas */

    /*
    ** do the initial iteration step
    */
    *iter =0; 
    strata[nused-1] =1;
    loglik[0] = coxfit6_iter(nvar, nused, method, beta, 1);
    loglik[1] = loglik[0];
    
    for(i = 0; i<nvar; i++){
      for(j = 0; j<nvar; j++){
        null_imat[i][j] = imat[i][j];
      }
    }

    /* am I done?
    **   update the betas and test for convergence
    */
    for (i=0; i<nvar; i++) /*use 'a' as a temp to save u0, for the score test*/
	a[i] = u[i];

    *flag= cholesky2(imat, nvar, toler);
    chsolve2(imat,nvar,a);        /* a replaced by  a *inverse(i) */

    temp=0;
    for (i=0; i<nvar; i++){
      temp +=  u[i]*a[i];
    }
    
    *sctest = temp;  /* score test */

    /*
    **  Never, never complain about convergence on the first step.  That way,
    **  if someone HAS to they can force one iter at a time.
    ** A non-finite loglik comes from exp overflow and requires almost
    **  malicious initial values.
    */
    if (maxiter==0 || isfinite(loglik[0])==0) {
	chinv2(imat,nvar);
	for (i=0; i<nvar; i++) {
	    beta[i] *= scale[i];  /*return to original scale */
	    u[i] /= scale[i];
	    imat[i][i] *= scale[i]*scale[i];
	    for (j=0; j<i; j++) {
		imat[j][i] *= scale[i]*scale[j];
		imat[i][j] = imat[j][i];
	    }
	}
	goto finish;
    }

    /*
    ** here is the main loop
    */
    loglik[1] = loglik[0];   /* loglik[1] contains the best so far */
    for (i=0; i<nvar; i++) {
	newbeta[i] = beta[i] + a[i];
    }

    halving = 0;       /* =1 when in the midst of "step halving" */
    for (*iter=1; *iter<= maxiter; (*iter)++) {
	R_CheckUserInterrupt();  
	newlk = coxfit6_iter(nvar, nused, method, newbeta, 0);

	/* am I done?
	**   test for convergence and then update beta
	*/
	*flag = cholesky2(imat, nvar, toler);

	notfinite = 0;
	for (i=0; i<nvar; i++) {
	    if (isfinite(u[i]) ==0) notfinite=2;     /* infinite score stat */
	    for (j=0; j<nvar; j++) {
		if (isfinite(imat[i][j]) ==0) notfinite =3; /*infinite imat */
	    }	
	}	
	if (isfinite(newlk) ==0) notfinite =4;
	
	if (notfinite==0 &&(fabs(1-(loglik[1]/newlk))<= eps)) { 
	    /* all done */
	    loglik[1] = newlk;
	    chinv2(imat, nvar);     /* invert the information matrix */

	    for (i=0; i<nvar; i++) {
		beta[i] = newbeta[i]*scale[i];
		u[i] /= scale[i];
		imat[i][i] *= scale[i]*scale[i];
		for (j=0; j<i; j++) {
		    imat[j][i] *= scale[i]*scale[j];
		    imat[i][j] = imat[j][i];
		}
	    }
	    if (halving) *flag= -2;
	    goto finish;
	}

	if (notfinite >0 || newlk < loglik[1])   {    
	    /*it is not converging ! */
	    halving++;  /* get more agressive when it doesn't work */
	    for (i=0; i<nvar; i++) 
		newbeta[i] = (newbeta[i] + halving*beta[i])/(halving + 1.0);
		    
	    }
	else {
	    halving=0;
	    loglik[1] = newlk;
	    chsolve2(imat,nvar,u);
	    for (i=0; i<nvar; i++) {
		beta[i] = newbeta[i];
		newbeta[i] = newbeta[i] +  u[i];
	    }	
	}
    }  /* return for another iteration */

    /*
    ** We end up here only if we ran out of iterations 
    **  recompute the last good version of imat and u,
    ** If maxiter =0 or 1, though, leave well enough alone.
    */
    if (maxiter > 1) 
	loglik[1] = coxfit6_iter(nvar, nused, method, beta, 0);
    chinv2(imat, nvar);
    for (i=0; i<nvar; i++) {
	beta[i] = beta[i]*scale[i];
	u[i] /= scale[i];
	imat[i][i] *= scale[i]*scale[i];
	for (j=0; j<i; j++) {
	    imat[j][i] *= scale[i]*scale[j];
	    imat[i][j] = imat[j][i];
	}
    }	
    *flag = 1000;
	
finish:
    /*
    ** create the output list
    */
    PROTECT(rlist= allocVector(VECSXP, 10));
    SET_VECTOR_ELT(rlist, 0, beta2);
    SET_VECTOR_ELT(rlist, 1, means2);
    SET_VECTOR_ELT(rlist, 2, u2);
    SET_VECTOR_ELT(rlist, 3, imat2);
    SET_VECTOR_ELT(rlist, 4, loglik2);
    SET_VECTOR_ELT(rlist, 5, sctest2);
    SET_VECTOR_ELT(rlist, 6, iter2);
    SET_VECTOR_ELT(rlist, 7, flag2);
    SET_VECTOR_ELT(rlist, 8, null_score2);
    SET_VECTOR_ELT(rlist, 9, null_imat2);
    

    /* add names to the objects */
    PROTECT(rlistnames = allocVector(STRSXP, 10));
    SET_STRING_ELT(rlistnames, 0, mkChar("coef"));
    SET_STRING_ELT(rlistnames, 1, mkChar("means"));
    SET_STRING_ELT(rlistnames, 2, mkChar("u"));
    SET_STRING_ELT(rlistnames, 3, mkChar("imat"));
    SET_STRING_ELT(rlistnames, 4, mkChar("loglik"));
    SET_STRING_ELT(rlistnames, 5, mkChar("sctest"));
    SET_STRING_ELT(rlistnames, 6, mkChar("iter"));
    SET_STRING_ELT(rlistnames, 7, mkChar("flag"));
    SET_STRING_ELT(rlistnames, 8, mkChar("null_score"));
    SET_STRING_ELT(rlistnames, 9, mkChar("null_imat"));
    setAttrib(rlist, R_NamesSymbol, rlistnames);

    unprotect(nprotect+2);
    return(rlist);
}

static double coxfit6_iter(int nvar, int nused,  int method, double *beta, int update_score) {
    int i, j, k, person;
    double  loglik =0;
    double  wtave;
    double  denom=0, zbeta, risk;
    double  temp2;
    int     ndead;  /* number of death obs at a time point */

    double  dtime;
    double  deadwt;  /*sum of case weights for the deaths*/
    double  denom2;  /* sum of weighted risk scores for the deaths*/
    int     nrisk;   /* number of subjects in the current risk set */
   
    for (i=0; i<nvar; i++) {
	u[i] =0;
	a2[i] =0;
	for(person=0; person < nused; person++){
	  if(update_score){
	    null_score[i][person] = 0;
	  }
	}
	for (j=0; j<nvar; j++) {
	    imat[i][j] =0 ;
	    cmat2[i][j] =0;
	    }
	}

    for (person=nused-1; person>=0; ) {
	if (strata[person] == 1) {
	    nrisk =0 ;  
	    denom = 0;
	    for (i=0; i<nvar; i++) {
		a[i] = 0;
		for (j=0; j<nvar; j++) cmat[i][j] = 0;
		}
	    }

	dtime = xtime[person];
	ndead =0; /*number of deaths at this time point */
	deadwt =0;  /* sum of weights for the deaths */
	denom2=0;  /* sum of weighted risks for the deaths */
  
  int start_person = person; /*id of first person in ties*/
  int end_person = person;
	while(person >=0 && xtime[person]==dtime) {
	    end_person = person;
	    /* walk through the this set of tied times */
	    nrisk++;
	    zbeta = offset[person];    /* form the term beta*z (vector mult) */
	    for (i=0; i<nvar; i++)
		zbeta += beta[i]*covar[i][person];
	    risk = exp(zbeta) * weights[person];
	    if (status[person] ==0) {
		denom += risk;
		/* a contains weighted sums of x, cmat sums of squares */
		for (i=0; i<nvar; i++) {
		    a[i] += risk*covar[i][person];
		    for (j=0; j<=i; j++)
			cmat[i][j] += risk*covar[i][person]*covar[j][person];
	        }
	    }	
	    else {
		ndead++;
		deadwt += weights[person];
		denom2 += risk;
		loglik += weights[person]*zbeta;

		for (i=0; i<nvar; i++) {
		    u[i] += weights[person]*covar[i][person];
		    a2[i] +=  risk*covar[i][person];
		    if(update_score){
		      null_score[i][person] += weights[person]*covar[i][person];
		    }
		    for (j=0; j<=i; j++)
			cmat2[i][j] += risk*covar[i][person]*covar[j][person];
		        }
	    }
	    person--;
	    if (person>=0 && strata[person]==1) break;  /*ties don't cross strata */
	    }

	if (ndead >0) {  /* we need to add to the main terms */
	    if (method==0 || ndead==1) { /* Breslow */
		denom += denom2;
		loglik -= deadwt* log(denom);
	   
		for (i=0; i<nvar; i++) {
		    a[i] += a2[i];
		    temp2= a[i]/ denom;  /* mean */
		    u[i] -=  deadwt* temp2;
		    if(update_score){
		      for(int iperson = start_person; iperson >= end_person; iperson--){
		        null_score[i][iperson] -= weights[iperson] * temp2;
		      }
		    }
		    for (j=0; j<=i; j++) {
			cmat[i][j] += cmat2[i][j];
			imat[j][i] += deadwt*(cmat[i][j] - temp2*a[j])/denom;
			imat[i][j] = imat[j][i];
			}
		    }
		}
	    else { /* Efron */
		/*
		** If there are 3 deaths we have 3 terms: in the first the
		**  three deaths are all in, in the second they are 2/3
		**  in the sums, and in the last 1/3 in the sum.  Let k go
		**  1 to ndead: we sequentially add a2/ndead and cmat2/ndead
		**  and efron_wt/ndead to the totals.
		*/
		wtave = deadwt/ndead;
		for (k=0; k<ndead; k++) {
		    denom += denom2/ndead;
		    loglik -= wtave* log(denom);
		    for (i=0; i<nvar; i++) {
			a[i] += a2[i]/ndead;
			temp2 = a[i]/denom;
			u[i] -= wtave *temp2;
			if(update_score){
			  for(int iperson = start_person; iperson >= end_person; iperson--){
			    null_score[i][iperson] -= weights[iperson]/ndead * temp2;
			  }
			}
			for (j=0; j<=i; j++) {
			    cmat[i][j] += cmat2[i][j]/ndead;
			    imat[j][i] += wtave*(cmat[i][j] - temp2*a[j])/denom;
			    imat[i][j] = imat[j][i];
			}	
		    }
		}
	    }	
	    for (i=0; i<nvar; i++) {
		a2[i]=0;
		for (j=0; j<nvar; j++) cmat2[i][j]=0;
	    }
	}
    }   /* end  of accumulation loop */
		
    /* A dummy line to stop a "set but never used" line in the compiler
    **  I don't use nrisk, but often used it in a print statement when debugging
    **  the code.   That may happen again 
    */
    nrisk = nrisk - 1;
    return(loglik);
}	
