/*
** Survival curves for a Cox model.  This routine counts up all
**  the totals that we need. The number at risk is a PITA in R code, but all
**  the rest of the compuations are simple there.  This creates more totals
**  than we currenly use, as an attempt to future proof the code.
** The routine uses a stacked input data set, one 'layer' per transition,
**  and outputs values for each transition, at each of the prespecified time
**  points. In a multistate model we need the values for all strata,
**  even if there was no event in some particular transtition at that time
**  point.  When there are extranal strata in the model, e.g.,
**  coxph(Surv(time, status) ~ age + strata(inst), lung), then this
**  routine will be called separatedly for individual strata, with otime the
**  event times in that stratum.  Ditto for a multistate model that had
**  external strata.
** 
**  otime:  vector of output times.  All the transitions will get reports at
**            these time points.  This fcn is called for all of the 
**            transitions at once, sorted by transition,
**            but called separately for any strata() groups.
**  y   :    survival response, two column
**  weight:  observation weight
**  sort2: sort index for the survival time
**  trans:   the data set is stacked: all the data for transition 1, then 
**            transition 2, etc for a multi-state model (th
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
**   6-7: censored: w1,w2
**   8-9: Efron number at risk
*/

#include <math.h>
#include "survS.h"
#include "survproto.h"
#include <stdio.h>

SEXP coxsurv1(SEXP otime2, SEXP y2,    SEXP weight2,  SEXP sort22, 
              SEXP trans2, SEXP xmat2, SEXP risk2) {
              
    int i, i2, k, person2, itrans;
    int nused, ntrans, ntime, irow, ii, jj;
    double *tstop, *status, *wt, *otime;
    int *trans;
    double dtime, meanwt;
    int  *sort2;
    double **xmat, *risk;  /* X matrix and risk score */
    int nvar;              /* number of covariates */
    double  *xsum1,  /* a weighted sum, for computing xbar */
	    *xsum2;
	    
    static const char *outnames[]={"ntrans", "count", 
				   "xbar", "xsum2", ""};
    SEXP rlist;
    double n[12];

    /* output variables */
    double  **rn, *rstrat;
    double  **rx1, **rx2;

    /* map the input data */
    otime = REAL(otime2);
    ntime = LENGTH(otime2);
    nused = nrows(y2);
    tstop = REAL(y2);
    status= tstop +nused;
    wt = REAL(weight2);
    sort2 = INTEGER(sort22);
    trans= INTEGER(trans2);
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
    
    /* Allocate memory for returned objects: ntime*ntrans copies of n, xsum1,
       and xsum2
    */
    PROTECT(rlist = mkNamed(VECSXP, outnames));
    irow = ntime*ntrans;
    rstrat = REAL(SET_VECTOR_ELT(rlist, 0, allocVector(REALSXP, 1)));
    rn = dmatrix(REAL(SET_VECTOR_ELT(rlist, 1, 
			    allocMatrix(REALSXP, irow, 10))), irow, 10);
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
    **  
    **  2. While tstop==otime, add them to the risk set, and count the 
    **   observation with repect to n3- n7.
    **
    ** In the code position2 is the current index into the sort2 vector,
    **   and i2 is the current value of sort2.
    */
    rstrat[0] = ntrans;   /* single element, number of transitions found */
    person2= nused-1;  
    irow = (ntime*ntrans);                /* row of output objects */
    for (ii =0; ii<ntrans; ii++) {
	itrans= trans[sort2[person2]];  /* current transition */
	for (k=0; k<10; k++) n[k] =0;
	for (k=0; k<nvar; k++) {xsum1[k] =0; xsum2[k] =0; }

	for (jj=ntime-1; jj>=0; jj--) { /* one by one through the times */
	    dtime = otime[jj];
	    for (k=3; k<8; k++) n[k]=0;  /* counts are only for this interval*/

	    /* Step 1 */
	    for(; person2 >=0 && trans[person2]==itrans; person2--) {
		i2 = sort2[person2];
		if (tstop[i2] <dtime) break;
		n[0]++;
		n[1] += wt[i2];
		n[2] += wt[i2] * risk[i2];
		for (k=0; k<nvar; k++) 
		    xsum1[k] += wt[i2]*risk[i2]*xmat[k][i2];

		if (status[i2]==0) {
		    /* count them as a 'censor' */
		    n[6]++;
		    n[7]+= wt[i2];
		}
		else if (tstop[i2]==dtime) {
		    /* step 2 */
printf("person=%d, i2=%d, dtime= %3.1f, n3=%3.1f\n", person2,i2, dtime, n[3]);
		    n[3]++;
		    n[4] += wt[i2];
		    n[5] += wt[i2]* risk[i2];
		    for (k=0; k<nvar; k++)
			xsum2[k] += wt[i2]*risk[i2]*xmat[k][i2];
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
	    for (k=0; k<10; k++) rn[k][irow] = n[k];
	    for (k=0; k<nvar; k++) {
		if (n[0]==0) rx1[k][irow]=  0;
		else         rx1[k][irow] = xsum1[k]/n[3];
		rx2[k][irow] = xsum2[k];
	    }
        } /* end of time points */

	/* walk past any data after the last selected time point */
	while(trans[person2]==itrans) person2--;
    }
    UNPROTECT(1);
    return(rlist);
}
    
