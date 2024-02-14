/*
** Compute the score residuals for a Cox model
**
** Input
**      y       matrix of time and status values
**      strata  = non-negative integer, unique value for each strata
**      covar2  the matrix of covariates, rows=variables, columns=subjects
**                (the S executive stores matrices in the Fortran ordering)
**      score   the vector of subject scores, i.e., exp(beta*z)
**      weights case weight
**      method  ==1 for efron method
**
** Output
**      resid   a matrix of the same shape as x
**
** Scratch
**      scratch,  from which a and a2 are carved
**
** Data must be sorted by strata, ascending time within strata.
**
** Updated 4/2023 to be O(np) instead of O(n^2 p): 
**   the score is sum (x_i -xbar(t)) (dN_i(t) - Y_i(t) risk_i lambda(t))
**   keep cumhaz = sum lambda(t) and xhaz = sum xbar(t) lambda(t) as running
**     totals, the second is a vector.  xbar requires loops that go from 
**     last time to first so that we can accumulate it.
** The main loop will go in "spurts": process all obs for this ending time,
**   all obs for the next ending time, etc. 
**   1, when we find a new time, for all obs at this time set initial resid to
**       risk * (x_i * cumhaz - xhaz)  (using the old cumhaz and xhaz)
**   2. find #events, hazard at this time point, mean at this time point, etc.
**     Update xbar, cumhaz and xhaz.
**   3. if obs is a death at this time, add (x_i - xbar(t))
**   4. at the end of a strata, subtract risk * (x_i*cumhaz - xhaz) for all
**     in the strata, and zero temporaries.
*/
#include <stdio.h>
#include "survS.h"
#include "survproto.h"

SEXP coxscore2(SEXP y2,     SEXP covar2,   SEXP strata2,
	       SEXP score2, SEXP weights2, SEXP method2) 
{
    int i,j, k, stratastart;
    int currentstrata;
    double temp;
    int n, nvar, method;
    double deaths, newtime;
    int dd;
    double *time, *status;
    double *score,  *weights;
    int *strata;
    double *a, *a2, *xhaz, xbar;
    double denom=0, e_denom;
    double risk;
    double **covar;
    double **resid;
    double hazard, cumhaz, meanwt;
    double downwt, temp2;
    SEXP   resid2;   /* the return matrix */

    n = nrows(y2);
    nvar = ncols(covar2);
    time = REAL(y2);
    status = time +n;
    strata = INTEGER(strata2);
    score  = REAL(score2);
    weights= REAL(weights2);
    method = asInteger(method2);

    /* scratch space */
    a = (double *) R_alloc(3*nvar, sizeof(double));
    a2 = a+nvar;
    xhaz = a2 + nvar;

    /*
    **  Set up the ragged arrays
    */
    covar=  dmatrix(REAL(covar2), n, nvar);
    PROTECT(resid2 = allocMatrix(REALSXP, n, nvar));
    resid=  dmatrix(REAL(resid2), n, nvar);

    denom=0; cumhaz=0;
    for (i=0; i<nvar; i++) {
	a2[i] =0;
	a[i] =0;
	xhaz[i] =0;
    }
    stratastart = n-1;
    currentstrata = strata[n-1];
    for (i=n-1; i >=0; ) {
	newtime = time[i];
	deaths =0; e_denom=0; meanwt =0;
	for (j=0; j< nvar; j++) a2[j] =0;
	
	for (; i>=0 && time[i]== newtime && strata[i] == currentstrata; i--) {
	    /* walk through any tied times */
	    risk = score[i] * weights[i];
	    denom += risk;
	    for (j=0; j<nvar; j++) {
		/* future accumulated risk that new entries don't get */
		resid[j][i] = score[i] * (covar[j][i]*cumhaz - xhaz[j]);
		a[j] += risk * covar[j][i]; /* running sum for covariates */
	    }
	    if (status[i]==1) {
		deaths++;
		e_denom += risk;
		meanwt += weights[i];
		for (j=0; j<nvar; j++) 
		    a2[j] += risk*covar[j][i];
	    }
	}
	
	if (deaths > 0) { /* update cumhaz and etc */
	    if (deaths <2 || method==0) {
		hazard = meanwt/denom;
		cumhaz += hazard;
		for (j=0; j<nvar; j++)  {
		    xbar = (a[j]/denom);     /* xbar for this variable */
		    xhaz[j] += xbar * hazard;
		    for (k=1+i; k<= i+ deaths; k++)
			resid[j][k] += covar[j][k] - xbar;
		}
	    }
	    else {  /* the harder case, Efron approx */
		/* If there are 3 deaths, the risk set includes all of
		**  them, then 2/3 of each, then 1/3 of each: think of it as
		**  3 separate deaths.  The censored people get all the cumhaz
		**  the deaths only a portion; we 'pre charge' them for the
		**  part of cumhaz and xhaz that they should not get at the
		**  end of the strata.
		*/
		meanwt /= deaths;
		for (dd=0; dd<deaths; dd++) {
		    downwt = dd/deaths;
		    temp = denom - downwt* e_denom;  /* working denominator */
		    hazard = meanwt/temp;
		    cumhaz += hazard;
		    for (j=0; j<nvar; j++) {
			xbar = (a[j] - downwt*a2[j])/ temp;
			xhaz[j] += xbar*hazard;
			for (k=1+i ; k<= i+ deaths; k++) {
			    temp2 = covar[j][k] - xbar;
			    resid[j][k] += temp2/deaths;
			    resid[j][k] += temp2 * score[k] * hazard * downwt;
			}
		    }
		}
	    }
	}
	
	if (i<0 || strata[i] != currentstrata) { /* end of a strata */
	    /* final work for each obs in the stratum */
	    for (k= stratastart; k> i; k--) {
		for (j=0; j < nvar; j++ )
		    resid[j][k] += score[k]* (xhaz[j] - covar[j][k]* cumhaz);
	    }
	    /* reset */
	    denom =0; cumhaz=0;
	    for (j=0; j<nvar; j++) {
		a[j] =0;
		xhaz[j] =0;
	    }
	    stratastart = i;
	    currentstrata= strata[i];
	}
    }

UNPROTECT(1);
return(resid2);
}
