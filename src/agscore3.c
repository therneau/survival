/*
** Do the score residuals
**
** Input
**      nx      number of subjects
**      nvarx   number of variables in the covariance matrix
**      y       matrix of start, stop, and event
**      strata  = unique non-negative for each stratum
**      covar2  the matrix of covariates, rows=variables, columns=subjects
**                (the S executive stores matrices in the Fortran ordering)
**      score   the vector of subject scores, i.e., exp(beta*z)
**      weights case weights
**      method  ==1 for efron approx
**      sort1 = sort order of y[,1], within strata
** Data is assumed to be sorted in y[,2] order, within strata. 
**  y[sort1[0],1] is the first entry time in strata 1, etc.
**
** Output
**      resid2  matrix of score residuals, same "shape" of matrix as covar2
**
** The score residual for subject i is \int (x- xbar(t)) dM_i(t), where
**  xbar is the weighted mean at time t and M is the martigale residual
**  process.  dM_i(t) = dN_i(t) - score[i] h(t)
**  h(t) is the increment in hazard = #events/ #at risk
**
** The key trick for efficiency is to only update the residual for any given
**  observation twice, once when they enter the risk set and again when they
**  leave, based on cumulative sums.  Go from largest to smallest time.
** The iteration goes in spurts: find a new time value, process all obs at that
**  time value, then go to the next.  
** 
** 1. Find the next time value, call it dtime.
** 2. For any who leave, Y[,1] >= dtime, 
**   a. finish the resid: subtract score* (x_i* cumhaz - xhaz) 
**   b. subtract them from denom and xsum   
** 3. If this is a new strata. 
**      a. For all in the prior strata who had not passed through step 2,
**          subtract score* (x_i* cumhaz - xhaz) from their residual
**      b. set cumhaz=0, xhaz= 0, xsum =0, denom =0
**      xhaz = cumulative sum of  hazard(t) * xbar(t), a vector
** 4. Over all observations with y[,2] = dtime
**    a. initialize their residual to  (x_i*cumhaz - xhaz)* score[i], before
**     updating cumhaz and xhaz
**    b. update xbar, cumhaz, xhaz with new information
**    c. for all deaths at the time, add (x_i- xbar) to their resid
**
** At the end of the loop finish up the final stratum.
**
*/
#include "survS.h"
#include "survproto.h"

SEXP agscore3(SEXP y2,     SEXP covar2,   SEXP strata2,
	      SEXP score2, SEXP weights2, SEXP method2, SEXP sort12) 
{
    int i,j, k;
    int n, nvar;
    int person, method;
    int currentstrata;
    double denom;
    double *a, *a2, *mean;
    int *strata;
    int i1, *sort1; 
    double *score, *weights;
    double e_denom;
    double risk, dtime;
    double hazard, meanwt;
    double  deaths, downwt;
    int dd;
    double *tstart, *tstop, *event;
    double **covar,
	   **resid;
    double cumhaz, *xhaz;
    double d2;
    double *mh1, *mh2, *mh3;
    SEXP resid2;  /* returned matrix */

    n = nrows(y2);
    nvar  = ncols(covar2);
    tstart = REAL(y2);
    tstop  = tstart +n;
    event = tstop + n;
    strata = INTEGER(strata2);
    score = REAL(score2);
    weights = REAL(weights2);
    method = asInteger(method2);
    sort1 = INTEGER(sort12);

    /* scratch space */
    a = (double *) R_alloc(7*nvar, sizeof(double));
    a2  = a+nvar;
    mean= a2 + nvar;
    mh1 = mean + nvar;
    mh2 = mh1 + nvar;
    mh3 = mh2 + nvar;
    xhaz = mh3 + nvar;

    /*
    **  Set up the ragged arrays
    */
    covar=  dmatrix(REAL(covar2), n, nvar);
    PROTECT(resid2 = allocMatrix(REALSXP, n, nvar));
    resid = dmatrix(REAL(resid2), n, nvar);
    for (i=0; i<n; i++) {
	for (k=0; k<nvar; k++) resid[k][i] =0.0;
    }	
    cumhaz =0; denom=0; 
    for (j=0; j<nvar; j++) {
	a[j] =0;
	xhaz[j] =0;
    }	
    i1 = n-1;
    currentstrata= strata[n-1];
    for (person=n-1; person>=0; ) {
	dtime = tstop[person];
	if (strata[person] != currentstrata) {
	    /* first obs of a new strata, finish off prior one */
	    for (; i1>=0 && sort1[i1] > person; i1--) {
		k = sort1[i1];
		for (j=0; j<nvar; j++)
		    resid[j][k] -= score[k]*(cumhaz*covar[j][k] - xhaz[j]);
	    }
	    /* rezero */
	    cumhaz =0; denom=0; 
	    for (j=0; j<nvar; j++) {
		a[j] =0;
		xhaz[j] =0;
	    }
	    currentstrata = strata[person];
	} else {
	    for(; i1>=0 && tstart[sort1[i1]] >= dtime; i1--) {
		k = sort1[i1];  /* observation to be removed from risk set */
		if (strata[k] != currentstrata) break;
		risk = score[k]*weights[k];
		denom -= risk;
		for (j=0; j<nvar; j++) {
		    resid[j][k] -= score[k]*(cumhaz*covar[j][k] - xhaz[j]);
		    a[j] -= risk *covar[j][k];
		}
	    }
	}
	
	/* count up over this time point */
	e_denom =0;
	meanwt =0;
	deaths =0;
	for (i=0; i<nvar; i++) a2[i]=0;

	for (; person >=0 && tstop[person] == dtime; person--) {
	    /* 
	    ** this next line is rare: the first obs of the next strata
	    **  has exactly the same time value as the last person
	    **  of the current strata
	    */
	    if (strata[person] != currentstrata) break;  
	    for (j=0; j<nvar; j++)
		resid[j][person] = (covar[j][person]*cumhaz - xhaz[j]) *
		    score[person];
	    risk = score[person] * weights[person];
	    denom += risk;             /* denominator of xbar(t) and hazard */
	    for (i=0; i<nvar; i++) {
		a[i] += risk*covar[i][person]; /* numerator of xbar(t) */
	    }

	    if (event[person]==1) {
		deaths++;
		e_denom += risk;
		meanwt += weights[person];
		for (i=0; i<nvar; i++)
		    a2[i] = a2[i] + risk*covar[i][person];
	    }
	}
	if (deaths >0) { /* update all the values */ 
	    if (deaths <2 || method==0) {
		/* easier case */
		hazard = meanwt/denom;
		cumhaz += hazard;
		for (i=0; i<nvar; i++) {
		    mean[i] = a[i]/denom;
		    xhaz[i] += mean[i]* hazard;
		    for (j=person+1; j<= person+ deaths; j++) {
			resid[i][j] += covar[i][j] - mean[i];
		    }
		}
	    }
	    else {
		/*
		** Efron case.  If there are k deaths, we treat it as k
		**  separate additions to the hazard.  For the second one
		**  each death has prob (k-1)/k of being present, then (k-2)/k,
		**  etc.  The idea is that the deaths actually occur in some
		**  order, we just don't know what that order is.
		** Say k=3 and h1, h2, h3 are the jumps in hazard.  The cumhaz
		**  and xhaz go up as we would expect.
		** The deaths get an addition  (x_i - xbar)/3 at each of the
		**  3 pseudo death times, and also a "look ahead" correction
		**  since they don't deserve the full increment of hazard.
		*/
		for (i=0; i<nvar; i++) {
		    mh1[i] =0;
		    mh2[i] =0;
		    mh3[i] =0;
		}
		meanwt /= deaths;  /* average weight of a death */
		for (dd=0; dd<deaths; dd++){
		    downwt = dd/deaths;
		    d2 = denom - downwt*e_denom;
		    hazard = meanwt/d2;
		    cumhaz += hazard;
		    for (i=0; i<nvar; i++) {
			mean[i] = (a[i] - downwt*a2[i])/ d2;
			xhaz[i] += mean[i]*hazard;
			mh1[i]  += hazard*downwt;
			mh2[i]  += mean[i] * hazard* downwt;
			mh3[i]  += mean[i]/deaths;
		    }
		}

		for (j=person+1; j<= person+ deaths; j++) {
		    for (i=0; i<nvar; i++) {
			resid[i][j] += (covar[i][j]-mh3[i]) + 
			    score[j]*(covar[i][j]*mh1[i] - mh2[i]);
		    }
		}
	    }	
	}
    }


    /* 
    ** finish those in the final stratum
    */
    for (; i1>=0; i1--) {
	k = sort1[i1];
	for (j=0; j<nvar; j++) 
	    resid[j][k] -= score[k] *(covar[j][k]*cumhaz  - xhaz[j]);
    }

    UNPROTECT(1);
    return(resid2);
}	
