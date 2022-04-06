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
#include "survS.h"
#include "survproto.h"

SEXP agscore2(SEXP y2,     SEXP covar2,   SEXP strata2,
	      SEXP score2, SEXP weights2, SEXP method2) 
{
    int i,k;
    int n, nvar;
    int person, method;
    double denom, time;
    double *a, *a2, *mean;
    int *strata;
    double *score, *weights;
    double e_denom;
    double risk;
    double hazard, meanwt;
    double  deaths, downwt;
    int dd;
    double *tstart, *tstop, *event;
    double **covar,
	   **resid;
    double temp1, temp2, d2;
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
   
    /* scratch space */
    a = (double *) R_alloc(6*nvar, sizeof(double));
    a2  = a+nvar;
    mean= a2 + nvar;
    mh1 = mean + nvar;
    mh2 = mh1 + nvar;
    mh3 = mh2 + nvar;

    /*
    **  Set up the ragged arrays
    */
    covar=  dmatrix(REAL(covar2), n, nvar);
    PROTECT(resid2 = allocMatrix(REALSXP, n, nvar));
    resid = dmatrix(REAL(resid2), n, nvar);
    for (i=0; i<n; i++) {
	for (k=0; k<nvar; k++) resid[k][i] =0.0;
    }	

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
	    time = tstop[person];
	    for (k=person; k<n; k++) {
		if (tstart[k] < time) {
		    risk = score[k] * weights[k];
		    denom += risk;
		    for (i=0; i<nvar; i++) {
			a[i] = a[i] + risk*covar[i][k];
		    }
		     if (tstop[k]==time && event[k]==1) {
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
	    if (deaths <2 || method==0) {
		/* easier case */
		hazard = meanwt/denom;
		for (i=0; i<nvar; i++) mean[i] = a[i]/denom;
		for (k=person; k<n; k++) {
		    if (tstart[k] < time) {
			risk = score[k];
			for (i=0; i<nvar; i++)
			    resid[i][k] -= (covar[i][k] -mean[i])*risk*hazard;
			if (tstop[k]==time) {
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
		    if (tstart[k] < time) {
			risk = score[k];
			if (tstop[k]==time && event[k]==1) {
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
		for ( ; tstop[person]==time; person++)
		    if (strata[person]==1) break;
	    }
	}
    }
UNPROTECT(1);
return(resid2);
}	
