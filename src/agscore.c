/*  $Id: agscore.c 11166 2008-11-24 22:10:34Z therneau $ */
/*
** Do the score residuals
**
** Input
**      nx      number of subjects
**      nvarx   number of variables in the covariance matrix
**      y       matrix of start, stop, and event
**      strata  =1 for the last obs of each strata
**      covar2  the matrix of covariates, rows=variables, columns=subjects
**                (the S executive stores matrices in the Fortran ordering)
**      score   the vector of subject scores, i.e., exp(beta*z)
**      weights case weights
**      method  ==1 for efron approx
**
** Output
**      resid2  matrix of score residuals, same "shape" of matrix as covar2
**
** Scratch
**      a       vector of length 3*nvar
*/
#include <stdio.h>
#include "survS.h"
#include "survproto.h"

void agscore(Sint   *nx,       Sint   *nvarx,      double *y,
	     double *covar2,   Sint   *strata,     double *score,
	     double *weights,  Sint   *method,     double *resid2, double *a)
    {
    int i,k;
    int n, nvar;
    int person;
    double denom, time;
    double *a2, *mean;
    double e_denom;
    double risk;
    double hazard, meanwt;
    double  deaths, downwt;
    int dd;
    double *start, *stop, *event;
    double **covar,
	   **resid;
    double temp1, temp2, d2;
    double *mh1, *mh2, *mh3;

    n = *nx;
    nvar  = *nvarx;
    start =y;
    stop  = y+n;
    event = y+(n+n);
    /*
    **  Set up the ragged arrays
    */
    covar=  dmatrix(covar2, n, nvar);
    resid = dmatrix(resid2, n, nvar);
    a2  = a+nvar;
    mean= a2 + nvar;
    mh1 = mean + nvar;
    mh2 = mh1 + nvar;
    mh3 = mh2 + nvar;

    for (person=0; person<n; ) {
	if (event[person]==0) person++;
	else {
	    /*
	    ** compute the mean over the risk set, also hazard at this time
	    */
	    denom =0;
	    e_denom =0;
	    meanwt =0;
	    deaths =0;
	    for (i=0; i<nvar; i++) {
		a[i] =0;
		a2[i]=0;
		}
	    time = stop[person];
	    for (k=person; k<n; k++) {
		if (start[k] < time) {
		    risk = score[k] * weights[k];
		    denom += risk;
		    for (i=0; i<nvar; i++) {
			a[i] = a[i] + risk*covar[i][k];
			}
		     if (stop[k]==time && event[k]==1) {
			deaths++;
			e_denom += risk;
			meanwt += weights[k];
			for (i=0; i<nvar; i++)
			    a2[i] = a2[i] + risk*covar[i][k];
			}
		     }
		if (strata[k]==1) break;
		}

	    /* add things in for everyone in the risk set*/
	    if (deaths <2 || *method==0) {
		/* easier case */
		hazard = meanwt/denom;
		for (i=0; i<nvar; i++) mean[i] = a[i]/denom;
		for (k=person; k<n; k++) {
		    if (start[k] < time) {
			risk = score[k];
			for (i=0; i<nvar; i++)
			    resid[i][k] -= (covar[i][k] -mean[i])*risk*hazard;
			if (stop[k]==time) {
			    person++;
			    if (event[k]==1)
				for (i=0; i<nvar; i++)
				    resid[i][k] += (covar[i][k] -mean[i]);
			    }
			}
		    if (strata[k]==1) break;
		    }
		}

	    else {
		/*
		** If there are 3 deaths, let m1, m2, m3 be the three
		**   weighted means,  h1, h2, h3 be the three hazard jumps.
		** Then temp1 = h1 + h2 + h3
		**      temp2 = h1 + (2/3)h2 + (1/3)h3
		**      mh1   = m1*h1 + m2*h2 + m3*h3
		**      mh2   = m1*h1 + (2/3)m2*h2 + (1/3)m3*h3
		**      mh3   = (1/3)*(m1+m2+m3)
		*/
		temp1=0;
		temp2=0;
		for (i=0; i<nvar; i++) {
		    mh1[i] =0;
		    mh2[i] =0;
		    mh3[i] =0;
		    }
		meanwt /= deaths;
		for (dd=0; dd<deaths; dd++){
		    downwt = dd/deaths;
		    d2 = denom - downwt*e_denom;
		    hazard = meanwt/d2;
		    temp1 += hazard;
		    temp2 += (1-downwt) * hazard;
		    for (i=0; i<nvar; i++) {
			mean[i] = (a[i] - downwt*a2[i])/ d2;
			mh1[i]  += mean[i] * hazard;
			mh2[i]  += mean[i] * (1-downwt) * hazard;
			mh3[i]  += mean[i]/deaths;
			}
		    }
		for (k=person; k<n; k++) {
		    if (start[k] < time) {
			risk = score[k];
			if (stop[k]==time && event[k]==1) {
			    for (i=0; i<nvar; i++) {
				resid[i][k] += covar[i][k] - mh3[i];
				resid[i][k] -= risk*covar[i][k]*temp2;
				resid[i][k] += risk* mh2[i];
				}
			    }
			else {
			    for (i=0; i<nvar; i++)
				resid[i][k] -= risk*(covar[i][k]*temp1 - mh1[i]);
			    }
			}
		    if (strata[k]==1) break;
		    }
		for ( ; stop[person]==time; person++)
		    if (strata[person]==1) break;
		}
	    }
	}
    }
