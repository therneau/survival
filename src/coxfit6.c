/*
** Cox regression fit, update to coxfit6, to use trust regions
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
**       step(nvar)            The proposed increment to beta
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

static double coxfit6_iter(int nvar, int nused, int method, double *beta);

SEXP coxfit6(SEXP maxiter2,  SEXP time2,   SEXP status2, 
	     SEXP covar2,    SEXP offset2, SEXP weights2,
	     SEXP strata2,   SEXP method2, SEXP eps2, 
	     SEXP toler2,    SEXP ibeta,    SEXP doscale2) {

    int i,j, person;
    double temp, temp2, temp3, temp4;
    double *newbeta, *scale;
    double halving =0, newlk;
    int notfinite;
    double *step, radius, egain, ratio;
    double *usave, *isave;

    /* copies of scalar input arguments */
    int     nused, nvar, maxiter;
    int     method;
    double  eps, toler;
    int *doscale;
   
    /* returned objects */
    SEXP imat2, means2, beta2, u2, loglik2;
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
    a = (double *) R_alloc(2*nvar*nvar + 6*nvar, sizeof(double));
    newbeta = a + nvar;
    a2 = newbeta + nvar;
    usave = a2 + nvar;
    scale = usave + nvar;
    step  = scale + nvar;
    cmat = dmatrix(step + nvar,   nvar, nvar);
    cmat2= dmatrix(step + nvar +nvar*nvar, nvar, nvar);
    isave= (double *) R_alloc(nvar*nvar, sizeof(double));

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

    /*
    ** Subtract the mean from each covar, as this makes the regression
    **  much more stable.
    ** The input variable doscale =0 for 0/1 indicator variables, they are
    **  left alone.
    ** Also set the intial step radius for the trust region, which will be
    **  min(3, 20/diff(range(x))).  The latter is to guarrantee that the
    **  the first beta is such that maxrisk/minrisk < exp(20).  The rationale
    **  for this choice is found in fail/trustregion.Rnw.
    */
    temp2=0;
    for (i=0; i<nused; i++) {
	temp2 += weights[i];
    }
    radius =3;  temp3=0; temp4=0;
    for (i=0; i<nvar; i++) {
	if (doscale[i]==0) {scale[i] = 1.0; means[i] = 0;}
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
		if (covar[i][person] < temp3) temp3 = covar[i][person];
		if (covar[i][person] > temp4) temp4 = covar[i][person];
		}
	    if (temp4 > temp3 &&  radius*(temp4-temp3) > 20)
		radius = 20/(temp4-temp3);  /* more conservative trust region*/
	    }
	}
 
    for (i=0; i<nvar; i++) beta[i] /= scale[i]; /*rescale initial betas */

    /*
    ** do the initial iteration step
    */
    *iter =0; 
    strata[nused-1] =1;
    loglik[0] = coxfit6_iter(nvar, nused, method, beta);
    loglik[1] = loglik[0];

    /* am I done?
    **   update the betas and test for convergence
    */
    for (i=0; i<nvar; i++) {
	step[i] = u[i];
	usave[i] = u[i];  /* needed if the first step goes badly awry */
    }
    
    *flag= cholesky2(imat, nvar, toler);
    chsolve2(imat,nvar, step);        /* u replaced by  u *inverse(i) */

    temp=0;
    for (i=0; i<nvar; i++)
	temp +=  u[i]*step[i];
    *sctest = temp;  /* score test */

    loglik[1] = loglik[0];   /* loglik[1] contains the best so far */
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
    ** Here is the main loop
    */
    for (*iter=1; *iter<= maxiter; (*iter)++) {
	R_CheckUserInterrupt();  

	if (*iter >1) {
	    /* compute the proposed step */
	    for (i=0; i<nvar; i++) step[i] = u[i];
	    chsolve2(imat, nvar, step);
	}
	temp =0;
	for (i=0; i<nvar; i++) temp += step[i]*step[i];
	temp = sqrt(temp);   /* length of proposed step */

	/*
	** Do the trust region checks.
	** step 1 is to check if the proposed step will lie within the
	**  trust region radius.  If so, use it.  If not, use a gradient
	**  step with length = radius (edge of trust region), and set halving=1
	**  as a flag.  
        ** Whichever is chosen, calculate the expected gain of the step.
	**  If the second order Taylor series that unerlies the Newton-Raphson
	**  iteration were perfect, the loglik will increase by s'U - s'Hs/2,
	**  where s is the step we are about to take and H= imat = second deriv
	*/
	egain =0;
	if (temp < radius) {
	    /* within the region */
	    halving =0;
	    for (i=0; i<nvar; i++) {
		newbeta[i] = beta[i] + step[i];
		egain += step[i]*u[i]/2;
	    }
	} else { /* constrained */
	    temp2 =0;
	    for (i=0; i<nvar; i++) temp2 += u[i]*u[i];
	    temp2 = sqrt(temp2);  /* length of the u vector */
	    for (i=0; i<nvar; i++) {
		step[i] = u[i]*radius/temp2;
		newbeta[i] = beta[i] + step[i];
	    }
	    for (i=0; i<nvar; i++) {
		egain += step[i]*u[i];
		temp3 = step[i];
		for (j=i+1; j<nvar; j++) {
		    /* we want (step %*% imat %*% step)/2, but remember
		       imat now contains the cholesky  */
		    temp3 += step[j]*imat[j][i];
		}
		egain -= temp3 * temp3 * imat[i][i]/2;
	    }
	    halving =1;  
	}
	    
	/* evaluate the loglik at this new point, and compute u and imat there*/
	newlk = coxfit6_iter(nvar, nused, method, newbeta);
	*flag = cholesky2(imat, nvar, toler);

	/* am I done?
	**   test for convergence and then update beta
	*/
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

	/*
	** Evaluate how well the step just taken worked, using the ratio of 
	**      (actual gain in loglik)/ (expected gain)
	** and use this to update the radius of the trust region.
        ** Then decide whether to accept the current increment or not
	*/
	ratio = (newlk - loglik[1])/egain;
	if (notfinite >0 || ratio < .25) radius = radius/4;
	else if (ratio > .75) radius = radius *2;

	if (notfinite ==0 && ratio > .1) {  /* accept */
	    for (i=0; i<nvar; i++) {
		beta[i] = newbeta[i];
		usave[i] = u[i];
	    }
	    for (i=0; i< nvar*nvar;  i++) 
		isave[i] = imat[0][i];
	    loglik[1] = newlk;
	} else {
	    /* reset and try again */
	    for (i=0; i< nvar; i++) u[i] = usave[i];
	    for (i=0; i< nvar*nvar; i++) imat[0][i] = isave[i];
	}
    }  /* return for another iteration */

    /*
    ** We end up here only if we ran out of iterations 
    **  recompute the last good version of imat and u,
    ** If maxiter =0 or 1, though, leave well enough alone.
    */
    if (maxiter > 1) 
	loglik[1] = coxfit6_iter(nvar, nused, method, beta);
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
    PROTECT(rlist= allocVector(VECSXP, 8));
    SET_VECTOR_ELT(rlist, 0, beta2);
    SET_VECTOR_ELT(rlist, 1, means2);
    SET_VECTOR_ELT(rlist, 2, u2);
    SET_VECTOR_ELT(rlist, 3, imat2);
    SET_VECTOR_ELT(rlist, 4, loglik2);
    SET_VECTOR_ELT(rlist, 5, sctest2);
    SET_VECTOR_ELT(rlist, 6, iter2);
    SET_VECTOR_ELT(rlist, 7, flag2);
    

    /* add names to the objects */
    PROTECT(rlistnames = allocVector(STRSXP, 8));
    SET_STRING_ELT(rlistnames, 0, mkChar("coef"));
    SET_STRING_ELT(rlistnames, 1, mkChar("means"));
    SET_STRING_ELT(rlistnames, 2, mkChar("u"));
    SET_STRING_ELT(rlistnames, 3, mkChar("imat"));
    SET_STRING_ELT(rlistnames, 4, mkChar("loglik"));
    SET_STRING_ELT(rlistnames, 5, mkChar("sctest"));
    SET_STRING_ELT(rlistnames, 6, mkChar("iter"));
    SET_STRING_ELT(rlistnames, 7, mkChar("flag"));
    setAttrib(rlist, R_NamesSymbol, rlistnames);

    unprotect(nprotect+2);
    return(rlist);
}

static double coxfit6_iter(int nvar, int nused,  int method, double *beta) {
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
	while(person >=0 && xtime[person]==dtime) {
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
		    for (j=0; j<=i; j++) {
			cmat[i][j] += cmat2[i][j];
			imat[j][i] += deadwt*(cmat[i][j] - temp2*a[j])/denom;
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
			for (j=0; j<=i; j++) {
			    cmat[i][j] += cmat2[i][j]/ndead;
			    imat[j][i] += wtave*(cmat[i][j] - temp2*a[j])/denom;
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
