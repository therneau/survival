/*
** Cox regression fit, replacement for coxfit2 in order
**   to be more frugal about memory: specificly that we 
**   don't make copies of the input data.
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
**       doscale      : 0=don't scale the X matrix, 1=scale the X matrix
**
**  returned parameters
**       means(nv)    : vector of column means of X
**       beta(nv)     :the vector of answers (at start contains initial est)
**       u(nv)        :score vector
**       imat(nv,nv)  :the variance matrix at beta=final
**                      (returned as a vector)
**       loglik(2)    :loglik at beta=initial values, at beta=final
**       sctest       :the score test at beta=initial
**       flag         :success flag  1000  did not converge
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

SEXP coxfit6(SEXP maxiter2,  SEXP time2,   SEXP status2, 
	     SEXP covar2,    SEXP offset2, SEXP weights2,
	     SEXP strata2,   SEXP method2, SEXP eps2, 
	     SEXP toler2,    SEXP ibeta,    SEXP doscale2) {
    int i,j,k, person;
    
    double **covar, **cmat, **imat;  /*ragged arrays */
    double  wtave;
    double *a, *newbeta;
    double *a2, **cmat2;
    double *scale;
    double  denom=0, zbeta, risk;
    double  temp, temp2;
    int     ndead;  /* actually, the sum of their weights */
    double  newlk=0;
    double  dtime, d2;
    double  deadwt;  /*sum of case weights for the deaths*/
    double  efronwt; /* sum of weighted risk scores for the deaths*/
    int     halving;    /*are we doing step halving at the moment? */
    int     nrisk;   /* number of subjects in the current risk set */
 
    /* copies of scalar input arguments */
    int     nused, nvar, maxiter;
    int     method;
    double  eps, toler;
    int doscale;

    /* vector inputs */
    double *time, *weights, *offset;
    int *status, *strata;
    
    /* returned objects */
    SEXP imat2, means2, beta2, u2, loglik2;
    double *beta, *u, *loglik, *means;
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
    doscale = asInteger(doscale2);

    time = REAL(time2);
    weights = REAL(weights2);
    offset= REAL(offset2);
    status = INTEGER(status2);
    strata = INTEGER(strata2);
    
    /*
    **  Set up the ragged arrays and scratch space
    **  Normally covar2 does not need to be duplicated, even though
    **  we are going to modify it, due to the way this routine was
    **  was called.  In this case NAMED(covar2) will =0
    */
    nprotect =0;
    if (NAMED(covar2)>0) {
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

    /*
    ** Subtract the mean from each covar, as this makes the regression
    **  much more stable.
    */
    for (i=0; i<nvar; i++) {
	temp=0;
	for (person=0; person<nused; person++) temp += covar[i][person];
	temp /= nused;
	means[i] = temp;
	for (person=0; person<nused; person++) covar[i][person] -=temp;
	if (doscale==1) {  /* and also scale it */
	    temp =0;
	    for (person=0; person<nused; person++) {
		temp += fabs(covar[i][person]);
	    }
	    if (temp > 0) temp = nused/temp;   /* scaling */
	    else temp=1.0; /* rare case of a constant covariate */
	    scale[i] = temp;
	    for (person=0; person<nused; person++)  covar[i][person] *= temp;
	    }
	}
    if (doscale==1) {
	for (i=0; i<nvar; i++) beta[i] /= scale[i]; /*rescale initial betas */
	}
    else {
	for (i=0; i<nvar; i++) scale[i] = 1.0;
	}

    /*
    ** do the initial iteration step
    */
    strata[nused-1] =1;
    loglik[1] =0;
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

	dtime = time[person];
	ndead =0; /*number of deaths at this time point */
	deadwt =0;  /* sum of weights for the deaths */
	efronwt=0;  /* sum of weighted risks for the deaths */
	while(person >=0 &&time[person]==dtime) {
	    /* walk through the this set of tied times */
	    nrisk++;
	    zbeta = offset[person];    /* form the term beta*z (vector mult) */
	    for (i=0; i<nvar; i++)
		zbeta += beta[i]*covar[i][person];
	    zbeta = coxsafe(zbeta);
	    risk = exp(zbeta) * weights[person];
	    denom += risk;

	    /* a is the vector of weighted sums of x, cmat sums of squares */
	    for (i=0; i<nvar; i++) {
		a[i] += risk*covar[i][person];
		for (j=0; j<=i; j++)
		    cmat[i][j] += risk*covar[i][person]*covar[j][person];
	        }

	    if (status[person]==1) {
		ndead++;
		deadwt += weights[person];
		efronwt += risk;
		loglik[1] += weights[person]*zbeta;

		for (i=0; i<nvar; i++) 
		    u[i] += weights[person]*covar[i][person];
		if (method==1) { /* Efron */
		    for (i=0; i<nvar; i++) {
			a2[i] +=  risk*covar[i][person];
			for (j=0; j<=i; j++)
			    cmat2[i][j] += risk*covar[i][person]*covar[j][person];
		        }
		    }
	        }
	    
	    person--;
	    if (strata[person]==1) break;  /*ties don't cross strata */
	    }


	if (ndead >0) {  /* we need to add to the main terms */
	    if (method==0) { /* Breslow */
		loglik[1] -= deadwt* log(denom);
	   
		for (i=0; i<nvar; i++) {
		    temp2= a[i]/ denom;  /* mean */
		    u[i] -=  deadwt* temp2;
		    for (j=0; j<=i; j++)
			imat[j][i] += deadwt*(cmat[i][j] - temp2*a[j])/denom;
		    }
		}
	    else { /* Efron */
		/*
		** If there are 3 deaths we have 3 terms: in the first the
		**  three deaths are all in, in the second they are 2/3
		**  in the sums, and in the last 1/3 in the sum.  Let k go
		**  from 0 to (ndead -1), then we will sequentially use
		**     denom - (k/ndead)*efronwt as the denominator
		**     a - (k/ndead)*a2 as the "a" term
		**     cmat - (k/ndead)*cmat2 as the "cmat" term
		**  and reprise the equations just above.
		*/
		for (k=0; k<ndead; k++) {
		    temp = (double)k/ ndead;
		    wtave = deadwt/ndead;
		    d2 = denom - temp*efronwt;
		    loglik[1] -= wtave* log(d2);
		    for (i=0; i<nvar; i++) {
			temp2 = (a[i] - temp*a2[i])/ d2;
			u[i] -= wtave *temp2;
			for (j=0; j<=i; j++)
			    imat[j][i] +=  (wtave/d2) *
				((cmat[i][j] - temp*cmat2[i][j]) -
					  temp2*(a[j]-temp*a2[j]));
		        }
		    }
		
		for (i=0; i<nvar; i++) {
		    a2[i]=0;
		    for (j=0; j<nvar; j++) cmat2[i][j]=0;
		    }
		}
	    }
	}   /* end  of accumulation loop */
    loglik[0] = loglik[1]; /* save the loglik for iter 0 */

    /* am I done?
    **   update the betas and test for convergence
    */
    for (i=0; i<nvar; i++) /*use 'a' as a temp to save u0, for the score test*/
	a[i] = u[i];

    *flag= cholesky2(imat, nvar, toler);
    chsolve2(imat,nvar,a);        /* a replaced by  a *inverse(i) */

    temp=0;
    for (i=0; i<nvar; i++)
	temp +=  u[i]*a[i];
    *sctest = temp;  /* score test */

    /*
    **  Never, never complain about convergence on the first step.  That way,
    **  if someone HAS to they can force one iter at a time.
    */
    for (i=0; i<nvar; i++) {
	newbeta[i] = beta[i] + a[i];
	}
    if (maxiter==0) {
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
    halving =0 ;             /* =1 when in the midst of "step halving" */
    for (*iter=1; *iter<= maxiter; (*iter)++) {
	newlk =0;
	for (i=0; i<nvar; i++) {
	    u[i] =0;
	    for (j=0; j<nvar; j++)
		imat[i][j] =0;
	    }

	/*
	** The data is sorted from smallest time to largest
	** Start at the largest time, accumulating the risk set 1 by 1
	*/
	for (person=nused-1; person>=0; ) {
	    if (strata[person] == 1) { /* rezero temps for each strata */
		denom = 0;
		nrisk =0;
		for (i=0; i<nvar; i++) {
		    a[i] = 0;
		    for (j=0; j<nvar; j++) cmat[i][j] = 0;
		    }
		}

	    dtime = time[person];
	    deadwt =0;
	    ndead =0;
	    efronwt =0;
	    while(person>=0 && time[person]==dtime) {
		nrisk++;
		zbeta = offset[person];
		for (i=0; i<nvar; i++)
		    zbeta += newbeta[i]*covar[i][person];
		zbeta = coxsafe(zbeta);
		risk = exp(zbeta) * weights[person];
		denom += risk;

		for (i=0; i<nvar; i++) {
		    a[i] += risk*covar[i][person];
		    for (j=0; j<=i; j++)
		    cmat[i][j] += risk*covar[i][person]*covar[j][person];
		    }

		if (status[person]==1) {
		    ndead++;
		    deadwt += weights[person];
		    newlk += weights[person] *zbeta;
		    for (i=0; i<nvar; i++) 
			u[i] += weights[person] *covar[i][person];
		    if (method==1) { /* Efron */
			efronwt += risk;
			for (i=0; i<nvar; i++) {
			    a2[i] +=  risk*covar[i][person];
			    for (j=0; j<=i; j++)
				cmat2[i][j] += risk*covar[i][person]*covar[j][person];
			    }   
		        }
	  	    }
		
		person--;
		if (strata[person]==1) break; /*tied times don't cross strata*/
	        }

	    if (ndead >0) {  /* add up terms*/
		if (method==0) { /* Breslow */
		    newlk -= deadwt* log(denom);
		    for (i=0; i<nvar; i++) {
			temp2= a[i]/ denom;  /* mean */
			u[i] -= deadwt* temp2;
			for (j=0; j<=i; j++)
			    imat[j][i] +=  (deadwt/denom)*
				(cmat[i][j] - temp2*a[j]);
		        }
    		    }
		else  { /* Efron */
		    for (k=0; k<ndead; k++) {
			temp = (double)k / ndead;
			wtave= deadwt/ ndead;
			d2= denom - temp* efronwt;
			newlk -= wtave* log(d2);
			for (i=0; i<nvar; i++) {
			    temp2 = (a[i] - temp*a2[i])/ d2;
			    u[i] -= wtave*temp2;
			    for (j=0; j<=i; j++)
				imat[j][i] +=  (wtave/d2)*
				    ((cmat[i][j] - temp*cmat2[i][j]) -
				    temp2*(a[j]-temp*a2[j]));
    		            }
    		        }

		    for (i=0; i<nvar; i++) { /*in anticipation */
			a2[i] =0;
			for (j=0; j<nvar; j++) cmat2[i][j] =0;
		        }
	            }
		}
	    }   /* end  of accumulation loop  */

	/* am I done?
	**   update the betas and test for convergence
	*/
	*flag = cholesky2(imat, nvar, toler);

	if (fabs(1-(loglik[1]/newlk))<= eps && halving==0) { /* all done */
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
	    goto finish;
	}

	if (*iter== maxiter) break;  /*skip the step halving calc*/

	if (newlk < loglik[1])   {    /*it is not converging ! */
		halving =1;
		for (i=0; i<nvar; i++)
		    newbeta[i] = (newbeta[i] + beta[i]) /2; /*half of old increment */
		}
	else {
	    halving=0;
	    loglik[1] = newlk;
	    chsolve2(imat,nvar,u);
	    j=0;
	    for (i=0; i<nvar; i++) {
		beta[i] = newbeta[i];
		newbeta[i] = newbeta[i] +  u[i];
	        }
	    }
	}   /* return for another iteration */

    /*
    ** We end up here only if we ran out of iterations 
    */
    loglik[1] = newlk;
    chinv2(imat, nvar);
    for (i=0; i<nvar; i++) {
	beta[i] = newbeta[i]*scale[i];
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
