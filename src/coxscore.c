/* $Id: coxscore.c 11357 2009-09-04 15:22:46Z therneau $
**
** Compute the score residuals for a Cox model
**
** Input
**      nx      number of subjects
**      nvarx   number of variables in the covariance matrix
**      y       matrix of time and status values
**      strata  =1 for the last obs of each strata
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
** Data must be sorted by strata, ascending time within strata, death before
**                      censor within time.
*/
#include <stdio.h>
#include "survS.h"
#include "survproto.h"

void coxscore(Sint   *nx,      Sint   *nvarx,    double *y, 
	      double *covar2,  Sint   *strata,   double *score, 
	      double *weights, Sint   *method,   double *resid2,
	      double *scratch)
    {
    int i,j, k;
    double temp;
    int n, nvar;
    double deaths;
    int dd;
    double *time, *status;
    double *a, *a2;
    double denom=0, e_denom;
    double risk;
    double **covar;
    double **resid;
    double hazard, meanwt;
    double downwt, temp2;
    double mean;

    n = *nx;
    nvar  = *nvarx;
    time = y;
    status = y+n;
    a = scratch;
    a2 = a+nvar;
    /*
    **  Set up the ragged array
    */
    covar=  dmatrix(covar2, n, nvar);
    resid=  dmatrix(resid2, n, nvar);

    e_denom=0;
    deaths=0;
    meanwt=0;
    for (i=0; i<nvar; i++) a2[i] =0;
    strata[n-1] =1;  /*failsafe */
    for (i=n-1; i >=0; i--) {
	if (strata[i]==1) {
	    denom =0;
	    for (j=0; j<nvar; j++) a[j] =0;
	    }

	risk = score[i] * weights[i];
	denom += risk;
	if (status[i]==1) {
	    deaths++;
	    e_denom += risk;
	    meanwt += weights[i];
	    for (j=0; j<nvar; j++) a2[j] += risk*covar[j][i];
	    }
	for (j=0; j<nvar; j++) {
	    a[j] += risk * covar[j][i];
	    resid[j][i] =0;
	    }

	if (deaths>0 && (i==0 || strata[i-1]==1 || time[i]!=time[i-1])){
	    /* last obs of a set of tied death times */
	    if (deaths <2 || *method==0) {
		hazard = meanwt/denom;
		for (j=0; j<nvar; j++)  {
		    temp = (a[j]/denom);     /* xbar */
		    for (k=i; k<n; k++) {
			temp2 = covar[j][k] - temp;
			if (time[k]==time[i] && status[k]==1)
				resid[j][k] += temp2;
			resid[j][k] -= temp2* score[k] * hazard;
			if (strata[k]==1) break;
			}
		    }
		}
	    else {  /* the harder case */
		meanwt /= deaths;
		for (dd=0; dd<deaths; dd++) {
		    downwt = dd/deaths;
		    temp = denom - downwt* e_denom;
		    hazard = meanwt/temp;
		    for (j=0; j<nvar; j++) {
			mean = (a[j] - downwt*a2[j])/ temp;
			for (k=i; k<n; k++) {
			    temp2 = covar[j][k] - mean;
			    if (time[k]==time[i] && status[k]==1) {
				resid[j][k] += temp2/deaths;
				resid[j][k] -= temp2 * score[k] * hazard *
						    (1 - downwt);
				}
			    else resid[j][k]-= temp2*score[k] * hazard;
			    if (strata[k]==1) break;
			    }
			}
		    }
		}
	    e_denom =0;
	    deaths =0;
	    meanwt =0;
	    for (j=0; j<nvar; j++)  a2[j] =0;
	    }
	}
    }
