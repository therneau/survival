/*
** Hazard curves for a Cox model.  This routine counts up all
**  the totals that we need (more than we need, actually).  The number at risk
**  is a PITA in R code.  Much of the rest of the computation can be done in R,
**  however.
**
**  y   :    survival response
**  weight:  observation weight
**  sort1, sort2: sort indices for the start and stop time
**  position: in a string of obs (1,2) (2,3) (3,4) (5,8) for a given subject,
**             1= start of a sequence, 2= end of sequence, 3= both
**            the above would be 1, 0, 2, 3
**  strata:  stratum for each obs
**  xmat2:   covariates
**  risk2:   risk score
**
**  output: a vector of unique death times, a vector of strata, a matrix of
**           counts (one row per time), and matrices of X sums, for all those
**           at risk and the current time.
**
**  For the weighted counts, number at risk != entries - exits.  Someone with 
**    a sequence of (1,2)(2,5)(5,6) will have 1 entry and 1 exit, but they might
**    have 3 changes of risk score due to time-dependent covariates.
**  So n1 has to count all the changes, while n2 and n3 (used for user printout)
**  only keep track of the first entry and final exit.
**
**  Let w1=1, w2= wt, w3= wt*risk.  The counts are
**   0-2: number at risk, w1, w2, w3
**   3-4: events: w1, w2
**   5-6: censor: w1, w2
**   7-9: censored endpoints: w1, w2, w3
**   10-11: entries: w1, w2
*/

#include <math.h>
#include "survS.h"
#include "survproto.h"
#include <stdio.h>

SEXP coxsurv1(SEXP y2, SEXP weight2,  SEXP sort12, SEXP sort22, 
              SEXP position2,   SEXP strata2, SEXP xmat2, SEXP risk2) {
              
    int i, i2, j2, k, person, person2, itime;
    int nused, istrat;
    int ntime;
    double *tstart=0, *stime, *status, *wt;
    int *strata;
    double dtime;
    int *sort1=0, *sort2;
    double **xmat, *risk;  /* X matrix and risk score */
    int nvar;              /* number of covariates */
    double  *xsum1,  /* a weighted sum, for computing xbar */
	    *xsum2;

    static const char *outnames[]={"time", "strata", "count", 
				   "xbar1", "xbar2", ""};
    SEXP rlist;
    double n[12];
    int *position;

    /* output variables */
    double  *rtime, **rn, *rstrat;
    double  **rx1, **rx2;

    /* map the input data */
    nused = nrows(y2);
    tstart = REAL(y2);
    stime = tstart + nused;
    status= stime +nused;
    wt = REAL(weight2);
    sort1 = INTEGER(sort12);
    sort2 = INTEGER(sort22);
    strata= INTEGER(strata2);
    position = INTEGER(position2);
    risk = REAL(risk2);
    nvar = ncols(xmat2);
    xmat = dmatrix(REAL(xmat2), nrows(xmat2), nvar);

    /* pass 1, get the number of unique event times, needed to alloc memory
    **  data is sorted by reverse time within strata
    */
    ntime =1; dtime=stime[sort2[0]]; istrat= strata[sort2[0]];
    for (i=1; i<nused; i++) {
	i2 = sort2[i];
	if (strata[i2] != istrat) {
	    ntime++;
	    istrat = strata[i2];
	}	
	else if (stime[i2] != dtime) ntime++;
    }  

    /* Allocate memory for the working matrices. */
    xsum1 = (double *) ALLOC(2*nvar, sizeof(double));
    xsum2 = xsum1 + nvar;
   
    /* Allocate memory for returned objects.  Essentially ntime copies of each
       of the above.  
       If ny=2 then the 'additions' of n1 are all zero, so no need to save them
    */
    PROTECT(rlist = mkNamed(VECSXP, outnames));
    rtime  = REAL(SET_VECTOR_ELT(rlist, 0, allocVector(REALSXP, ntime)));
    rstrat = REAL(SET_VECTOR_ELT(rlist, 1, allocVector(REALSXP, ntime)));
    rn = dmatrix(REAL(SET_VECTOR_ELT(rlist, 2, 
			    allocMatrix(REALSXP, ntime, 12))), ntime, 12);
    rx1 = dmatrix(REAL(SET_VECTOR_ELT(rlist, 3, 
			      allocMatrix(REALSXP, ntime, nvar))), ntime, nvar);
    rx2 = dmatrix(REAL(SET_VECTOR_ELT(rlist, 4, 
			      allocMatrix(REALSXP, ntime, nvar))), ntime, nvar);
						     
    R_CheckUserInterrupt();  /*check for control-C */

    /* now add up all the sums */
    person = 0; person2=0 ;   /* person2 tracks start times */
    istrat = strata[sort2[0]];  /* initial stratum */
    for (itime=ntime-1; itime>=0; itime--) {
	i2 = sort2[person];
	if (person==0 || strata[i2] != istrat) {
	    if (person>0) {
		/* catch up on the entries */
		j2 = sort2[person2];
		for (; tstart[j2] >= dtime && strata[j2]==istrat; person2++) {
		    n[10]++;
		    n[11] += wt[j2];
		    j2 = sort1[person2];
		}
printf("itime=%d, n= %4.2f %4.2f\n", itime, n[10], n[11]);
		rn[10][itime+1] = n[10];
		rn[11][itime+1] = n[11];
	    }
            /* new stratum, zero everything */
	    for (k=0; k<12; k++) n[k] =0;
	    for (k=0; k < nvar; k++) xsum1[k] =0;
	    istrat = strata[i2];
	}
	dtime = stime[i2];    /*current time point */
	rtime[itime] = dtime;
	rstrat[itime] = istrat;
	for (k=3; k<12; k++) n[k]=0;  /* zero all but cumulative n */

	while(person< nused && stime[i2]==dtime && strata[i2]==istrat) {
	    for (k=0; k<nvar; k++) xsum2[k] =0;

	    n[0]++;
	    n[1] += wt[i2];
	    n[2] += wt[i2] * risk[i2];
	    for (k=0; k<nvar; k++) xsum1[k] += wt[i2]*risk[i2]*xmat[k][person];

	    if (status[i2] > 0) { /* an event */
		for (k=0; i<nvar; k++) 
		    xsum2[k] += wt[i2]*risk[i2]* xmat[k][person];
		n[3]++;
		n[4]+= wt[i2];
		if (position[i2] >1) {
		    n[7]++;
		    n[8]+= wt[i2];
		    n[9] += wt[i2]*risk[i2];
		}
	    }
	    
	    if (position[i2] > 1) {  /* 2 or 3 = end of a subject's series */
		n[5]++;
		n[6] += wt[i2];
	    }
	    person++;
	    i2 = sort2[person];
	}				

	j2 = sort1[person2];
	while(person2 < nused && tstart[j2] >= dtime && strata[j2]==istrat) {
	    /* remove subjects */
	    n[0]--;
	    if (n[0] ==0) {
		n[1] =0;
		n[2] =0;
		for (k=0; k<nvar; k++) xsum1[k] =0;
	    }
	    else {
		n[1] -= wt[j2];
		n[2] -= wt[j2]*risk[j2];
		for (k=0; k<nvar; k++) xsum1[k] -=xmat[k][j2] * wt[j2]* risk[j2];
	    }
		    
	    /* count entries */
	    if (position[j2]==1 || position[j2]==3) {
		n[10]++;
		n[11] += wt[j2];
	    }
		
	    person2++;
	    j2 = sort1[person2];	
	}
/*    printf("itime=%d, istrat=%d, dtime=%f, n=%3.1f %3.1f %3.1f %3.1f %3.1f \n",
	       itime, istrat, dtime, n[0], n[3], n[5], n[7], n[10]);
*/	
	/* copy results to the R structures */
	for (k=0; k<12; k++) rn[k][itime] = n[k];
	for (k=0; k<nvar; k++){
	    rx1[k][itime] = xsum1[k]/n[3];
	    rx2[k][itime] = xsum2[k]/n[3];
	}
    }

    j2 = sort2[person2];
    for (; person2 < nused; person2++) {
	n[10]++;
	n[11] += wt[j2];
    }
    rn[10][0] = n[10];
    rn[11][0] = n[11];

    UNPROTECT(1);
    return(rlist);
}
    
