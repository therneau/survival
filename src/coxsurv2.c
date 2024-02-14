/*
** Survival curves for a Cox model.  This routine counts up all
**  the totals that we need. The number at risk is a PITA in R code, but all
**  the rest of the compuations are simple there.  This creates more totals
**  than we currenly use, as an attempt to future proof the code.
** 
** coxsurv1: (time,status) for multistate
** coxsurv2: (time1, time2, status) for multistate
** coxsurv3: (time, status) for single state
** coxsurv4: (time1, time2, status) for single state
**
**  otime:  vector of output times.  All the transitions will get reports at
**            these time points.  (To be useful, the times should be a superset
**            of the event times, since the parent routine will do cumsums.)
**  y   :    survival response
**  weight:  observation weight
**  sort1, sort2: sort indices for the start and stop time
**  sindex: in a sequence of obs (1,2) (2,3) (3,4) (5,8) for a given subject,
**             'sindex' would be 1,    0,    2,    3
**             1= start of a sequence, 2= end of sequence, 3= both
**  trans:   the data set is stacked: all the data for transition 1, then 
**            transition 2, etc for a multi-state model
**  xmat2:   covariates
**  risk2:   risk score
**
**  output: A matrix of counts, with one row per time, per stratum, and
**           matrices with xbar for those at risk, and the sum of x for
**	     terminal events at the current time.
**
**  For the weighted counts, number at risk != entries - exits.  Someone with 
**    a sequence of (1,2)(2,5)(5,6) will have 1 entry and 1 exit, but they might
**    have 3 changes of risk score due to time-dependent covariates.
**  n0-3 has to count all the changes, while n8-n9 (only used in printout)
**    keep track of the final exit, and 3-7 and 10-11 refer only to a given
**    timepoint.
**
**  Let w1=1, w2= wt, w3= wt*risk.  The counts n[] are
**   0-2: number at risk: w1, w2, w3
**   3-5: events: w1, w2, w3
**   6-7: censored endpoints: w1,w2
**   8-9: Efron sums 1 and 2
**  10-11: censored counts w1 and w2
*/

#include <math.h>
#include "survS.h"
#include "survproto.h"
#include <stdio.h>

SEXP coxsurv2(SEXP otime2, SEXP y2, SEXP weight2,  SEXP sort12, SEXP sort22, 
              SEXP sindex2,  SEXP trans2, SEXP xmat2, SEXP risk2) {
              
    int i, i1, i2, k, person1, person2, itrans;
    int nused, ntrans, ntime, irow, ii, jj;
    double *tstart=0, *tstop, *status, *wt, *otime;
    int *trans;
    double dtime, meanwt;
    int *sort1=0, *sort2;
    double **xmat, *risk;  /* X matrix and risk score */
    int nvar;              /* number of covariates */
    double  *xsum1,  /* a weighted sum, for computing xbar */
	    *xsum2;
	    
    int *atrisk;

    static const char *outnames[]={"ntrans", "count", 
				   "xbar", "xsum2", ""};
    SEXP rlist;
    double n[12];
    int *sindex;

    /* output variables */
    double  **rn, *rstrat;
    double  **rx1, **rx2;

    /* map the input data */
    otime = REAL(otime2);
    ntime = LENGTH(otime2);
    nused = nrows(y2);
    tstart = REAL(y2);
    tstop = tstart + nused;
    status= tstop +nused;
    wt = REAL(weight2);
    sort1 = INTEGER(sort12);
    sort2 = INTEGER(sort22);
    trans= INTEGER(trans2);
    sindex = INTEGER(sindex2);
    risk = REAL(risk2);
    nvar = ncols(xmat2);
    xmat = dmatrix(REAL(xmat2), nrows(xmat2), nvar);

    /* pass 1, count the number of transitions, needed to alloc memory
    **  data is sorted by time within trans
    */
    itrans= trans[0];
    ntrans=1;
    for (i=1; i<nused; i++) {
	i2 = sort2[i];
	if (trans[i2] != itrans) {
	    ntrans++;
	    itrans = trans[i2];
	}	
    }  
 
    /* Allocate memory for the working matrices. */
    xsum1 = (double *) ALLOC(2*nvar, sizeof(double));
    xsum2 = xsum1 + nvar;
    atrisk = (int *) ALLOC(nused, sizeof(int));
    for (i=0; i<nused; i++) atrisk[i] =0;
    
    /* Allocate memory for returned objects: ntime*ntrans copies of n, xsum1,
       and xsum2
    */
    PROTECT(rlist = mkNamed(VECSXP, outnames));
    irow = ntime*ntrans;
    rstrat = REAL(SET_VECTOR_ELT(rlist, 0, allocVector(REALSXP, 1)));
    rn = dmatrix(REAL(SET_VECTOR_ELT(rlist, 1, 
			    allocMatrix(REALSXP, irow, 12))), irow, 12);
    rx1 = dmatrix(REAL(SET_VECTOR_ELT(rlist, 2, 
			      allocMatrix(REALSXP, irow, nvar))), irow, nvar);
    rx2 = dmatrix(REAL(SET_VECTOR_ELT(rlist, 3, 
			      allocMatrix(REALSXP, irow, nvar))), irow, nvar);
						     
    R_CheckUserInterrupt();  /*check for control-C */

    /* now add up all the sums 
    **  All this is done backwards in time.  The logic is a bit easier, and
    **   the computation is numerically more stable (fewer subtractions).
    **  One by one for the desired output times "otime".
    **  1. While tstop > otime, walk backwards adding subjects to the risk
    **     set when we cross their ending time.  Also add them to the "censored"
    **     count.  Don't add any who will be removed before otime to either
    **     count, however.
    **    Set atrisk=1 if they are added.   
    **  
    **  2. While tstop==otime, add them to the risk set, and count the 
    **   observation with repect to n3- n7.
    **
    **  3. While tstart > otime, remove any obs currently at risk from n0-n2.
    **
    ** In the code position1, position2 are current indices into the sort1
    **   and sort2 vectors, and i1/i2 are the current values of those vectors.
    */
    rstrat[0] = ntrans;   /* single element, number of transitions found */
    person1 = nused-1; person2= nused-1;   /* person1 tracks start times */    
    irow = (ntime*ntrans);                /* row of output objects */
    for (ii =0; ii<ntrans; ii++) {
	itrans= trans[sort2[person2]];  /* current transition */
	for (k=0; k<3; k++) n[k] =0;
	for (k=0; k<nvar; k++) {xsum1[k] =0; xsum2[k] =0; }

	for (jj=ntime-1; jj>=0; jj--) { /* one by one through the times */
	    dtime = otime[jj];
	    for (k=3; k<12; k++) n[k]=0;  /* counts are only for this interval*/

	    /* Step 1 */
	    for(; person2 >=0 && trans[person2]==itrans; person2--) {
		i2 = sort2[person2];
		if (tstop[i2] <dtime) break;
		if (tstart[i2] < dtime) {
		    /* add them to the risk set */
		    /* if otime were (10, 20) an interval (14,15) will never
		       be added  */
		    atrisk[i2] =1;
		    n[0]++;
		    n[1] += wt[i2];
		    n[2] += wt[i2] * risk[i2];
		    for (k=0; k<nvar; k++) 
			xsum1[k] += wt[i2]*risk[i2]*xmat[k][i2];

		    if (sindex[i2]>1 && status[i2]==0) {
			/* count them as a 'censor' */
			n[10]++;
			n[11]+= wt[i2];
		    }
		}

		if (tstop[i2]==dtime && status[i2]>0) {
		    /* step 2  (tstart < tstop BTW)*/
		    n[3]++;
		    n[4] += wt[i2];
		    n[5] += wt[i2]* risk[i2];
		    for (k=0; k<nvar; k++)
			xsum2[k] += wt[i2]*risk[i2]*xmat[k][i2];
		    if (sindex[i2] >1) {
			n[6]++;
			n[7] += wt[i2];
		    }
		}
	    }

	    /* Step 3 */
	    for (; person1 >=0 && trans[person1]==itrans; person1--) {
		i1 = sort1[person1];
		if (tstart[i1] < dtime) break;
		if (atrisk[i1]) {  /* remove them from risk set */
		    n[0]--; 
		    if (n[0] ==0) {
			n[1] =0;
			n[2] =0;
			for (k=0; k<nvar; k++) xsum1[k] =0;
		    }
		    else {
			n[1] -= wt[i1];
			n[2] -= wt[i1]*risk[i1];
			for (k=0; k<nvar; k++) 
			    xsum1[k] -= xmat[k][i1] * wt[i1]* risk[i1];
		    }
		}
	    }

	    /* Compute the Efron number at risk */
	    if (n[3] <=1) {   /* only one event */
		n[8]= n[2];
		n[9] = n[2]*n[2];
	    }		
	    else {
		meanwt = n[5]/(n[3]*n[3]);  /* average weight of deaths /n */
		for (k=0; k<n[3]; k++) {
		    n[8] += n[2] - k*meanwt;
		    n[9] += (n[2] -k*meanwt)*(n[2] - k*meanwt);
		}
		n[8] /= n[3];
		n[9] /= n[3];
	    }		

	    /* save the results */
	    irow--;
	    if (irow <0) Rprintf("irow error in coxsurv2.  This should never happen: please contact package author\n");
	    for (k=0; k<12; k++) rn[k][irow] = n[k];
	    for (k=0; k<nvar; k++) {
		if (n[0]==0) rx1[k][irow]=  0;
		else         rx1[k][irow] = xsum1[k]/n[3];
		rx2[k][irow] = xsum2[k];
	    }
        } /* end of time points */

	/* walk past any data after the last selected time point */
	while(trans[person2]==itrans) person2--;
	while(trans[person1]==itrans) person1--;
    }
    UNPROTECT(1);
    return(rlist);
}
    
