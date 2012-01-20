/*
** Cox regression fit, replacement for coxfit2 in order
**   to be more frugal about memory: specificly that the 
**   x and y matrices are not copied back to the output
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
**         actually, I package loglik onward into one object
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
	     SEXP covar2,    SEXP offset2, SEXP weights,
	     SEXP strata2,   SEXP method2, SEXP eps2, 
	     SEXP toler2,    SEXP ibeta,    SEXP doscale2) {
    int i,j,k, person;
    int     iter;
    
    double **covar, **cmat, **imat;  /*ragged arrays */
    double *mark, *wtave;
    double *a, *newbeta;
    double *a2, **cmat2;
    double  denom=0, zbeta, risk;
    double  temp, temp2;
    double  ndead;
    double  newlk=0;
    double  d2, efron_wt;
    int     halving;    /*are we doing step halving at the moment? */
    int     nrisk;   /* number of subjects in the current risk set */
 
    /* copies of scalar input arguments */
    int     nused, nvar, maxiter;
    int     method;
    double  eps, toler;
    int doscale;
    
    /* returned objects */
    SEXP imat2, means2, beta2, u2, loglik2;
    double *beta, *u, *loglik;
    SEXP rlist, rlistnames;
    int nprotect;  /* number of protect calls I have issued */

    nused = LENGTH(offset2);
    nvar  = ncol(covar2);
    method = asInteger(method2);
    maxiter = asInteger(maxiter2);
    eps  = asReal(eps2);     /* convergence criteria */
    toler = asReal(toler2);  /* tolerance for cholesky */
    doscale = asInteger(doscale2);
    
    /*
    **  Set up the ragged arrays and scratch space
    **  We need to duplicate covar, since we will be scaling it
    */
    if (doscale) covar= dmatrix(REAL(duplicate(covar2)), nused, nvar);
    else  covar =  dmatrix(REAL(covar2), nused, nvar);
    PROTECT(imat2 = allocVector(REALSXP, nvar*nvar)); 
    imat = dmatrix(REAL(imat2),  nvar, nvar);
    a = R_alloc(2*nvar*nvar + 2*nused + 3*nvar);
    newbeta = a + nvar;
    a2 = newbeta + nvar;
    mark = a2 + nvar;
    wtave= mark + nused;
    cmat = dmatrix(wtave + nvar,   nvar, nvar);
    cmat2= dmatrix(wtave + nvar +nvar*nvar, nvar, nvar);

    /* 
    ** create output variables
    */ 
    PROTECT(beta2 = duplicate(ibeta));
    beta = REAL(beta2);
    PROTECT(means2 = allocVector(REALSXP, nvar));
    means = REAL(means2);
    PROTECT(u2 = allocVector(REALSXP, nvar));
    u = REAL(u2);
    PROTECT(loglik2 = allocVector(REALSXP, 5)); /* loglik, sctest, flag,maxiter*/
    loglik = REAL(loglik2);
    nprotect = 5;

    /*
    **   Mark(i) contains the number of tied deaths at this point,
    **    for the first person of several tied times. It is zero for
    **    the second and etc of a group of tied times.
    **   Wtave contains the average weight for the deaths
    */
    temp=0;
    j=0;
    for (i=nused-1; i>0; i--) {
	if ((time[i]==time[i-1]) & (strata[i-1] != 1)) {
	    j += status[i];
	    temp += status[i]* weights[i];
	    mark[i]=0;
	    }
	else  {
	    mark[i] = j + status[i];
	    if (mark[i] >0) wtave[i]= (temp+ status[i]*weights[i])/ mark[i];
	    temp=0; j=0;
	    }
	}
    mark[0]  = j + status[0];
    if (mark[0]>0) wtave[0] = (temp +status[0]*weights[0])/ mark[0];

    /*
    ** Subtract the mean from each covar, as this makes the regression
    **  much more stable.
    */
    for (i=0; i<nvar; i++) {
	for (person=0; person<nused; person++) temp += covar[i][person];
	temp /= nused;
	means[i] = temp;
	if (doscale) {
	    for (person=0; person<nused; person++) covar[i][person] -=temp;
	    temp =0;
	    for (person=0; person<nused; person++) temp += abs(covar[i][person]);
	    temp = nused/temp;   /* scaling */
	    scale[i] = temp;
	    for (person=0; person<nused; person++)  covar[i][person]*temp;
	    }
	}
    if (doscale) {
	for (i=0; i<nvar; i++) beta[i] /= scale;
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
	for (j=0; j<nvar; j++)
	    imat[i][j] =0 ;
	}

    efron_wt =0;
    for (person=nused-1; person>=0; person--) {
	if (strata[person] == 1) {
	    nrisk =0 ;  
	    denom = 0;
	    for (i=0; i<nvar; i++) {
		a[i] = 0;
		a2[i]=0 ;
		for (j=0; j<nvar; j++) {
		    cmat[i][j] = 0;
		    cmat2[i][j]= 0;
		    }
		}
	    }
	
	nrisk++;
	zbeta = offset[person];    /* form the term beta*z   (vector mult) */
	for (i=0; i<nvar; i++)
	    zbeta += beta[i]*covar[i][person];
	zbeta = coxsafe(zbeta);
	risk = exp(zbeta) * weights[person];
	if (method==2) r2[nrisk] = risk;
	denom += risk;
	efron_wt += status[person] * risk;  /*sum(denom) for tied deaths*/

	for (i=0; i<nvar; i++) {
	    a[i] += risk*covar[i][person];
	    for (j=0; j<=i; j++)
		cmat[i][j] += risk*covar[i][person]*covar[j][person];
	    }

	if (status[person]==1) {
	    loglik[1] += weights[person]*zbeta;
	    for (i=0; i<nvar; i++) {
		    u[i] += weights[person]*covar[i][person];
		a2[i] +=  risk*covar[i][person];
		for (j=0; j<=i; j++)
		    cmat2[i][j] += risk*covar[i][person]*covar[j][person];
		}
	    }
	if (mark[person] >0) {  /* once per unique death time */
	    ndead = mark[person];
	    if (method==0 || ndead==1) { /* Breslow */
		loglik[1] -= ndead * wtave[person]* log(d2);
		for (i=0; i<nvar; i++) {
		    temp2= a[i]/ denom;  /* mean */
		    u[i] -= wtave*ndead * temp2;
		    for (j=0; j<=i; j++)
			imat[j][i] +=  ndead *wtave[person]*(
				 (cmat[i][j] - temp*cmat2[i][j]) /denom -
					  temp2*a[j]/denom);
		    }
		}
	    else { /* Efron */
		for (k=0; k<ndead; k++) {
		    temp = (double)k * method / ndead;
		    d2= denom - temp*efron_wt;
		    loglik[1] -= wtave[person] * log(d2);
		    for (i=0; i<nvar; i++) {
			temp2 = (a[i] - temp*a2[i])/ d2;
			u[i] -= wtave[person] *temp2;
			for (j=0; j<=i; j++)
			    imat[j][i] +=  wtave[person]*(
				 (cmat[i][j] - temp*cmat2[i][j]) /d2 -
					  temp2*(a[j]-temp*a2[j])/d2);
		        }
		    }
		}

	    efron_wt =0;
	    for (i=0; i<nvar; i++) {
		a2[i]=0;
		for (j=0; j<nvar; j++)  cmat2[i][j]=0;
		}
	    }
	}   /* end  of accumulation loop */

    loglik[0] = loglik[1];   /* save the loglik for iteration zero  */

    /* am I done?
    **   update the betas and test for convergence
    */
    for (i=0; i<nvar; i++) /*use 'a' as a temp to save u0, for the score test*/
	a[i] = u[i];

    *flag= cholesky2(imat, nvar, *tol_chol);
    chsolve2(imat,nvar,a);        /* a replaced by  a *inverse(i) */

    temp=0;
    for (i=0; i<nvar; i++)
	temp +=  u[i]*a[i];
    loglik[3] = temp;  /* score test */
    /*
    **  Never, never complain about convergence on the first step.  That way,
    **  if someone HAS to they can force one iter at a time.
    */
    for (i=0; i<nvar; i++) {
	newbeta[i] = beta[i] + a[i];
	}
    if (*maxiter==0) {
	chinv2(imat,nvar);
	for (i=1; i<nvar; i++) {
	    beta[i] *= scale;  /*return to original scale */
	    imat[i][i] *= scale[i]*scale[i];
	    for (j=0; j<i; j++) {
		imat[j][i] *= scale[i]*scale[j];
		imat[i][j] = imat[j][i];
		}
	return;  
	}

    /*
    ** here is the main loop
    */
    halving =0 ;             /* =1 when in the midst of "step halving" */
    for (iter=1; iter<=*maxiter; iter++) {
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
	for (person=nused-1; person>=0; person--) {
	    if (strata[person] == 1) { /* rezero temps for each strata */
		efron_wt =0;
		denom = 0;
		nrisk =0;
		for (i=0; i<nvar; i++) {
		    a[i] = 0;
		    a2[i]=0 ;
		    for (j=0; j<nvar; j++) {
			cmat[i][j] = 0;
			cmat2[i][j]= 0;
			}
		    }
		}

	    nrisk++;
	    zbeta = offset[person];
	    for (i=0; i<nvar; i++)
		zbeta += newbeta[i]*covar[i][person];
	    zbeta = coxsafe(zbeta);
	    risk = exp(zbeta ) * weights[person];
	    denom += risk;
	    efron_wt += status[person] * risk;  /* sum(denom) for tied deaths*/

	    for (i=0; i<nvar; i++) {
		a[i] += risk*covar[i][person];
		for (j=0; j<=i; j++)
		    cmat[i][j] += risk*covar[i][person]*covar[j][person];
		}

	    if (status[person]==1) {
		newlk += weights[person] *zbeta;
		for (i=0; i<nvar; i++) {
		    u[i] += weights[person] *covar[i][person];
		    a2[i] +=  risk*covar[i][person];
		    for (j=0; j<=i; j++)
			cmat2[i][j] += risk*covar[i][person]*covar[j][person];
		    }
		}

	    if (mark[person] >0) {  /* once per unique death time */
		ndead = mark[person];
		if (method==0 || ndead==1) { /* Breslow */
		    loglik[1] -= ndead * wtave[person]* log(d2);
		    for (i=0; i<nvar; i++) {
			temp2= a[i]/ denom;  /* mean */
			u[i] -= wtave*ndead * temp2;
			for (j=0; j<=i; j++)
			    imat[j][i] +=  ndead *wtave[person]*(
				(cmat[i][j] - temp*cmat2[i][j]) /denom -
				temp2*a[j]/denom);
    		        }
    		    }
		else if (method==1) { /* Efron */
		    for (k=0; k<ndead; k++) {
			temp = (double)k * method / ndead;
			d2= denom - temp*efron_wt;
			loglik[1] -= wtave[person] * log(d2);
			for (i=0; i<nvar; i++) {
			    temp2 = (a[i] - temp*a2[i])/ d2;
			    u[i] -= wtave[person] *temp2;
			    for (j=0; j<=i; j++)
				imat[j][i] +=  wtave[person]*(
				    (cmat[i][j] - temp*cmat2[i][j]) /d2 -
				    temp2*(a[j]-temp*a2[j])/d2);
    		            }
    		        }
    		    }
		else {  /* Exact calculation */	
		    }

		efron_wt =0;
		for (i=0; i<nvar; i++) {
		    a2[i]=0;
		    for (j=0; j<nvar; j++)  cmat2[i][j]=0;
		    }
		}
	    }   /* end  of accumulation loop  */

	/* am I done?
	**   update the betas and test for convergence
	*/
	*flag = cholesky2(imat, nvar, *tol_chol);

	if (fabs(1-(loglik[1]/newlk))<=*eps && halving==0) { /* all done */
	    loglik[1] = newlk;
	    chinv2(imat, nvar);     /* invert the information matrix */
	    for (i=1; i<nvar; i++) {
		beta[i] = newbeta[i]*scale[i];
		imat[i][i] *= scale[i]*scale[i];
		for (j=0; j<i; j++) {
		    imat[j][i] *= scale[i]*scale[j];
		    imat[i][j] = imat[j][i];
		    }
		loglik[2] = iter;
	    goto finish;
	    }

	if (iter==*maxiter) break;  /*skip the step halving calc*/

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

    loglik[1] = newlk;
    chinv2(imat, nvar);
    for (i=1; i<nvar; i++) {
	beta[i] = newbeta[i]*scale[i];
	imat[i][i] *= scale[i]*scale[i];
	for (j=0; j<i; j++) {
	    imat[j][i] *= scale[i]*scale[j];
	    imat[i][j] = imat[j][i];
	    }
	}
    loglik[4] = 1000;

finish:
    /*
    ** create the output list
    */
    PROTECT(rlist= allocVector(VECSXP, 5));
    SET_VECTOR_ELT(rlist, 0, beta2);
    SET_VECTOR_ELT(rlist, 1, means2);
    SET_VECTOR_ELT(rlist, 2, u2);
    SET_VECTOR_ELT(rlist, 3, imat2);
    SET_VECTOR_ELT(rlist, 4, loglik2);

    /* add names to the objects */
    PROTECT(rlistnames = allocVector(STRSXP, 5));
    SET_STRING_ELT(rlistnames, 0, mkChar("coef"));
    SET_STRING_ELT(rlistnames, 1, mkChar("means"));
    SET_STRING_ELT(rlistnames, 2, mkChar("u"));
    SET_STRING_ELT(rlistnames, 3, mkChar("imat"));
    SET_STRING_ELT(rlistnames, 4, mkChar("loglik"));
    setAttrib(rlist, R_NamesSymbol, rlistnames);

    unprotect(nprotect+2);
    return(rlist);
    }
