/*
** Hazard curves for a Cox model.  This routine counts up all
**  the totals that we need (more than we need, actually).  The number at risk
**  is a PITA in R code.  Much of the rest of the computation can be done in R,
**  however.
**
** coxsurv1: (time,status) for multistate
** coxsurv2: (time1, time2, status) for multistate
** coxsurv3: (time, status) for single state
** coxsurv4: (time1, time2, status) for single state
**
**  y   :    survival response
**  weight:  observation weight
**  sort2: sort indices for the stop time
**  strata:  stratum for each obs
**  xmat2:   covariates
**  risk2:   risk score
**
**  output: a vector of unique death times, a vector of strata, a matrix of
**           counts (one row per time), and matrices of X sums, for all those
**           at risk and the current time.
**
**  Let w1=1, w2= wt, w3= wt*risk.  The counts are
**   0-2: number at risk, w1, w2, w3
**   3-4: events: w1, w2
**   5-6: censor: w1, w2
*/

#include <math.h>
#include "survS.h"
#include "survproto.h"
#include <stdio.h>

SEXP coxsurv3(SEXP y2, SEXP weight2,  SEXP sort22, 
              SEXP strata2, SEXP xmat2, SEXP risk2) {
              
    int i, i2, k, person, itime;
    int nused, istrat;
    int ntime;
    double *stime, *status, *wt;
    int *strata;
    double dtime;
    int *sort2;
    double **xmat, *risk;  /* X matrix and risk score */
    int nvar;              /* number of covariates */
    double  *xsum1,  /* a weighted sum, for computing xbar */
	    *xsum2;

    static const char *outnames[]={"time", "strata", "count", 
				   "xbar1", "xbar2", ""};
    SEXP rlist;
    double n[8];

    /* output variables */
    double  *rtime, **rn, *rstrat;
    double  **rx1, **rx2;

    /* map the input data */
    nused = nrows(y2);
    stime = REAL(y2);
    status= stime +nused;
    wt = REAL(weight2);
    sort2 = INTEGER(sort22);
    strata= INTEGER(strata2);
    risk = REAL(risk2);
    nvar = ncols(xmat2);
    xmat = dmatrix(REAL(xmat2), nrows(xmat2), nvar);

     /* pass 1, get the number of unique times, needed to alloc memory
    **  data is sorted by time within strata. The same time in 2 strata
    **  counts as 2 separate times.
    */
    ntime =1; dtime=stime[sort2[0]]; istrat= strata[sort2[0]];
    for (i=1; i<nused; i++) {
	i2 = sort2[i];
	if (strata[i2] != istrat) {
	    ntime++;
	    istrat = strata[i2];
	}	
	else if (stime[i2] != dtime) ntime++;
	dtime = stime[i2];
    }  

    /* Allocate memory for the working matrices. */
    xsum1 = (double *) ALLOC(2*nvar, sizeof(double));
    xsum2 = xsum1 + nvar;
    /* Allocate memory for returned objects.  Essentially ntime copies of each
       of n, xsum1, xsum2.
    */
    PROTECT(rlist = mkNamed(VECSXP, outnames));
    rtime  = REAL(SET_VECTOR_ELT(rlist, 0, allocVector(REALSXP, ntime)));
    rstrat = REAL(SET_VECTOR_ELT(rlist, 1, allocVector(REALSXP, ntime)));
    rn = dmatrix(REAL(SET_VECTOR_ELT(rlist, 2, 
			    allocMatrix(REALSXP, ntime, 7))), ntime, 7);
    rx1 = dmatrix(REAL(SET_VECTOR_ELT(rlist, 3, 
			      allocMatrix(REALSXP, ntime, nvar))), ntime, nvar);
    rx2 = dmatrix(REAL(SET_VECTOR_ELT(rlist, 4, 
			      allocMatrix(REALSXP, ntime, nvar))), ntime, nvar);
						     
    R_CheckUserInterrupt();  /*check for control-C */

    for (k=0; k<7; k++) n[k] =0;
    for (k=0; k < nvar; k++) xsum1[k] =0;
    /* now add up all the sums */
    person = nused -1;
    istrat = strata[sort2[person]];  /* initial stratum */

    for (itime=ntime-1; itime>=0; itime--) {
	i2 = sort2[person];
	if (strata[i2] != istrat) {
            /* new stratum, zero everything */
	    for (k=0; k<7; k++) n[k] =0;
	    for (k=0; k < nvar; k++) xsum1[k] =0;
	    istrat = strata[i2];
	}
	dtime = stime[i2];    /*current time point */
	rtime[itime] = dtime;
	rstrat[itime] = istrat;
	for (k=3; k<8; k++) n[k]=0;  /* zero all but cumulative n */

	while(person>0 && stime[i2]==dtime && strata[i2]==istrat) {
	    for (k=0; k<nvar; k++) xsum2[k] =0;
	    n[0]++;
	    n[1] += wt[i2];
	    n[2] += wt[i2] * risk[i2];
	    for (k=0; k<nvar; k++) xsum1[k] += wt[i2]*risk[i2]*xmat[k][person];

	    if (status[i2] > 0) { /* an event */
		for (k=0; k<nvar; k++) 
		    xsum2[k] += wt[i2]*risk[i2]* xmat[k][person];
		n[3]++;
		n[4]+= wt[i2];
	    } else {
		n[5]++;
		n[6]+= wt[i2];
	    }
	    
	    person--;
	    i2 = sort2[person];
	}				

	/* copy results to the R structures */
	for (k=0; k<nvar; k++){
	    rx1[k][itime] = xsum1[k]/n[3];
	    rx2[k][itime] = xsum2[k]/n[3];
	}
    }

    UNPROTECT(1);
    return(rlist);
}
    
