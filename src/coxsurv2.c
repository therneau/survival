/*
** Survival curves for a Cox model.  This routine counts up all
**  the totals that we need (more than we need, actually).  The number at risk
**  is a PITA in R code.  Much of the rest of the computation can be done in R,
**  however.
**
**  otime:  vector of output times.  All the strata will get reports at these
**            time points.  (Normally this is called for all of the states at
**            once, each state being a stratum.)
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
**  output: A matrix of counts, with one row per time, per stratum, and
**           matrices with xbar for those at risk, and the sum of x for
**	     terminal events at the current time.
**
**  For the weighted counts, number at risk != entries - exits.  Someone with 
**    a sequence of (1,2)(2,5)(5,6) will have 1 entry and 1 exit, but they might
**    have 3 changes of risk score due to time-dependent covariates.
**  So n0-3 has to count all the changes, while n8-n11 (used only for
**   printout) only keep track of the first entry and final exit, and 3-7
**   are events at a given timepoint.
**
**  Let w1=1, w2= wt, w3= wt*risk.  The counts are
**   0-2: number at risk, w1, w2, w3
**   3-4: events: w1, w2
**   5-7: censored endpoints: w1, w2, w3
**   8-9: censor: w1, w2
**  10-11: entry: w1, w2
*/

#include <math.h>
#include "survS.h"
#include "survproto.h"
#include <stdio.h>

SEXP coxsurv2(SEXP otime2, SEXP y2, SEXP weight2,  SEXP sort12, SEXP sort22, 
              SEXP position2,  SEXP strata2, SEXP xmat2, SEXP risk2) {
              
    int i, i2, j2, k, person, person2, istrat;
    int nused, nstrat, ntime, irow, ii, jj;
    double *tstart=0, *tstop, *status, *wt, *otime;
    int *strata;
    double dtime;
    int *sort1=0, *sort2;
    double **xmat, *risk;  /* X matrix and risk score */
    int nvar;              /* number of covariates */
    double  *xsum1,  /* a weighted sum, for computing xbar */
	    *xsum2;
    int *atrisk;

    static const char *outnames[]={"nstrat", "count", 
				   "xbar", "xsum2", ""};
    SEXP rlist;
    double n[12];
    int *position;

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
    strata= INTEGER(strata2);
    position = INTEGER(position2);
    risk = REAL(risk2);
    nvar = ncols(xmat2);
    xmat = dmatrix(REAL(xmat2), nrows(xmat2), nvar);

    /* pass 1, count the number of strata, needed to alloc memory
    **  data is sorted by time within strata
    */
    istrat= strata[sort2[0]];
    nstrat=1;
    for (i=1; i<nused; i++) {
	i2 = sort2[i];
	if (strata[i2] != istrat) {
	    nstrat++;
	    istrat = strata[i2];
	}	
    }  
 
    /* Allocate memory for the working matrices. */
    xsum1 = (double *) ALLOC(2*nvar, sizeof(double));
    xsum2 = xsum1 + nvar;
    atrisk = (int *) ALLOC(nused, sizeof(int));
    for (i=0; i<nused; i++) atrisk[i] =0;
    
    /* Allocate memory for returned objects: ntime*nstrat copies of n, xsum1,
       and xsum2
    */
    PROTECT(rlist = mkNamed(VECSXP, outnames));
    irow = ntime*nstrat;
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
    **     It that was their first entry, add them to an entry count for the
    **     prior otime (n9, n10).
    **
    ** At the end of a strata add any remaining subject to n10-11, which 
    **   works on a lag
    */
    rstrat[0] = nstrat;
    person = nused-1; person2= nused-1;   /* person2 tracks start times */
    irow = (ntime*nstrat); n[11]=0; n[10]=0;
    for (ii =0; ii<nstrat; ii++) {
	istrat= strata[sort2[person]];  /* current stratum */
	for (k=0; k<3; k++) n[k] =0;
	for (k=0; k<nvar; k++) {xsum1[k] =0; xsum2[k] =0; }

	for (jj=ntime-1; jj>=0; jj--) { /* one by one through the times */
	    dtime = otime[jj];
	    for (k=3; k<10; k++) n[k]=0;  /* counts are only for this interval*/

	    /* Step 1 */
	    i2 = sort2[person];
	    while(tstop[i2] >= dtime && person>=0 && strata[person]==istrat) {
		if (tstart[i2] <= dtime) {
		    /* add them to the risk set */
		    atrisk[i2] =1;
		    n[0]++;
		    n[1] += wt[i2];
		    n[2] += wt[i2] * risk[i2];
		    for (k=0; k<nvar; k++) 
			xsum1[k] += wt[i2]*risk[i2]*xmat[k][person];
		}

		if (position[i2]>1) {
		    /* count them as an 'censor' */
		    n[8]++;
		    n[9]+= wt[j2];
		}
		
		if (tstop[i2]==dtime && status[i2]>0) {
		    /* step 2 */
		    n[3]++;
		    n[4] += wt[i2];
		    if (position[i2] >1) {
			n[5]++;
			n[6] += wt[i2];
			n[7] += wt[i2]*risk[i2];
		    }
		}
		person--;
		i2 = sort2[person];
	    }
	    /*printf("A ii=%d, jj=%d, person=%d, n0=%3.1f, n3=%3.1f, n8=%3.1f\n", 
		   ii, jj,  person, n[0], n[3], n[8]);*/

	    /* Step 3 */
	    j2 = sort1[person2];
	    while(tstart[j2] >= dtime && person2 >=0 && strata[person2]==istrat){
		if (atrisk[j2]) {  /* remove them from risk set */
		    n[0]--; 
		    if (n[0] ==0) {
			n[1] =0;
			n[2] =0;
			for (k=0; k<nvar; k++) xsum1[k] =0;
		    }
		    else {
			n[1] -= wt[j2];
			n[2] -= wt[j2]*risk[j2];
			for (k=0; k<nvar; k++) 
			    xsum1[k] -=xmat[k][j2] * wt[j2]* risk[j2];
		    }
		}

		/* additions that are "in the gap" between two otimes */
		if (position[j2] ==1 || position[j2]==3) { /* count as a entry */
		    n[10]++;
		    n[11]+= wt[j2];
		}
		/*printf("person2=%d, j2=%d, position=%d, n10=%3.1f\n",
		       person2, j2, position[j2], n[10]); */
		person2--;
		j2 = sort2[person2];
	    }
	    /*printf(" B person=%d, n0=%3.1f, n10=%3.1f\n", person, n[0], n[10]);*/

	    /* save the results */
	    if (jj< (ntime-1)) {
		/* save the entries as part of the later otime */
		rn[10][irow] = n[10];
		rn[11][irow] = n[11];
		n[10]=0;
		n[11]=0;
	    }   

	    irow--;
	    for (k=0; k<10; k++) rn[k][irow] = n[k];
	    for (k=0; k<nvar; k++) {
		if (n[0]==0) rx1[k][irow]=0;
		else         rx1[k][irow] = xsum1[k]/n[3];
		rx2[k][irow] = xsum2[k];
	    }
        } /* end of time points */

 	while(person2>=0 && strata[person2]==istrat) {
	    /* censors after the last time are listed at the last time */
	    n[10]++;
	    n[11]+= wt[j2];
	    person2--;
	    j2 = sort2[j2];
	    rn[10][irow] = n[10];
	    rn[11][irow] = n[11];
	}	
	
    }
    UNPROTECT(1);
    return(rlist);
}
    
