/*
** Hazard curves for a Cox model.  This routine counts up all
**  the totals that we need (more than we need, actually).  The number at risk
**  is a PITA in R code.  Much of the rest of the computation could be done in 
**  R, however.
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
**  risk2:   risk score = exp(linear.predictor) from the Cox fit
**  efron2:  1 if the Efron approx, 0 for Breslow
**
**  output: a vector of unique death times, a vector of strata, a matrix of
**           counts (one row per time), 
**           the matrix of mean X values (ntime x nvar)
**           and the score residuals (n x nvar)
**
**  Let w1=1, w2= wt, w3= wt*risk.  The counts are
**   0-2: number at risk; w1, w2, w3 weighting
**   3-5: events: w1, w2, w3
**   6  : modified n[2] for the Efron approx
**
**  The score residual calculations are lifted from coxscore2.c, see that
**   for comments about the algorithm.
*/

#include <math.h>
#include "survS.h"
#include "survproto.h"
#include <stdio.h>

SEXP coxsurv3(SEXP y2, SEXP xmat2,  SEXP strata2, 
              SEXP risk2, SEXP weight2, SEXP sort22,
	      SEXP efron2) {
              
    int i, i2, j, k, itime, dd;
    int nused, istrat, stratastart;
    int efron, ntime;
    double *stime, *status, *wt;
    int *strata;
    double dtime, meanwt;
    double hazard, cumhaz;
    double temp, temp2, downwt, tmean;
    int *sort2;
    double **xmat, *risk;  /* X matrix and risk score */
    int nvar;              /* number of covariates */
    double  *xsum1;  /* a weighted sum, for computing xbar */
    double  *xsum2;  /* needed for the Efron approx */
    double  *xhaz, *xmean;  /* last is a temp needed for Efron approx */
    static const char *outnames[]={"time", "strata", "count", 
				   "xbar", "sresid", ""};
    SEXP rlist;
    double n[7];

    /* output variables */
    double  *rtime, **rn;
    int     *rstrat;
    double  **xbar, **resid;

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
    efron = asInteger(efron2);

    /* pass 1, get the number of unique event times, needed to alloc memory
    **  data is sorted by time within strata. The same event time in 2 strata
    **  counts as 2 separate times.
    */
    ntime =0; dtime=stime[sort2[0]] -1; istrat= strata[sort2[0]];
    for (i=0; i<nused; i++) {
	i2 = sort2[i];
	if (strata[i2] != istrat) {
	    istrat = strata[i2];  /* new stratum */
	    dtime  = stime[i2] -1; /* smallest in new strata -1 */
	}
	
	if (status[i2] == 1 && stime[i2] != dtime) {
	    ntime++;
	    dtime = stime[i2];
	}
    }  

    /* Allocate memory for the working matrices and returned objects */
    xsum1 = (double *) ALLOC(4*nvar, sizeof(double));
    xsum2 = xsum1 + nvar;
    xhaz  = xsum2 + nvar;
    xmean = xhaz  + nvar;

    PROTECT(rlist = mkNamed(VECSXP, outnames));
    rtime  = REAL(SET_VECTOR_ELT(rlist, 0, allocVector(REALSXP, ntime)));
    rstrat = INTEGER(SET_VECTOR_ELT(rlist, 1, allocVector(INTSXP, ntime)));
    rn = dmatrix(REAL(SET_VECTOR_ELT(rlist, 2, 
			    allocMatrix(REALSXP, ntime, 7))), ntime, 7);
    xbar = dmatrix(REAL(SET_VECTOR_ELT(rlist, 3, 
			      allocMatrix(REALSXP, ntime, nvar))), ntime, nvar);
    resid = dmatrix(REAL(SET_VECTOR_ELT(rlist, 4, 
			      allocMatrix(REALSXP, nused, nvar))), nused, nvar);
						     
    R_CheckUserInterrupt();  /*check for control-C */

    for (j=0; j<7; j++) n[j] =0;
    for (j=0; j < nvar; j++) {
	xsum1[j] =0;
	xhaz[j]  =0;
    }
    cumhaz=0;
    /* now add up all the sums */
    itime = ntime -1;
    stratastart = nused -1;  /* start point for the first stratum */
    istrat = strata[sort2[nused-1]];  /* initial stratum */
    for (i= nused -1; i >=0; ) {
	i2 = sort2[i];
	dtime = stime[i2];
        for (j=0; j< nvar; j++) xsum2[j] =0;
	for (j=3; j< 7; j++) n[j]=0;

	if (strata[i2] != istrat) {
	    /* finish score residuals for prior stratum */
	    for (k= stratastart; k> i; k--) {
                for (j=0; j < nvar; j++ )
                    resid[j][k] += risk[k]* (xhaz[j] - xmat[j][k]* cumhaz);
            }

            /* new stratum, zero totals */
	    for (j=0; j<3; j++) n[j] =0;
	    for (j=0; j < nvar; j++) {
		xsum1[j] =0;
		xhaz[j]  =0; 
	    }
	    cumhaz=0;
	    istrat = strata[i2];
	    stratastart = i;
	}
	
       for (; i>=0; i--) {
	   i2 = sort2[i];
	   if (stime[i2]!= dtime || strata[i2] != istrat) break;
            /* walk through any tied times */
 	    n[0]++;
	    n[1] += wt[i2];
	    n[2] += wt[i2] * risk[i2];
	    for (j=0; j<nvar; j++) {
		xsum1[j] += wt[i2]*risk[i2]* xmat[j][i2];
		resid[j][i2] = risk[i2]*(xmat[j][i2]*cumhaz - xhaz[j]);
	    }

	    if (status[i2] ==1) {
		n[3]++;
		n[4]+= wt[i2];
		n[5] += wt[i2]* risk[i2];
		for (j=0; j<nvar; j++)
                    xsum2[j] += risk[i2] * wt[i2] *xmat[j][i2];
 	    } 
       }

       if (n[3] >0) { /* if any deaths */
	   /* update terms for score residuals */    
	   if (n[3] <2 || efron ==0) {
	       hazard = n[4]/n[2];
	       cumhaz += hazard;
	       for (j=0; j<nvar; j++)  {
		   xmean[j] = xsum1[j]/n[2];
		   xhaz[j] += xmean[j] * hazard;
		   for (k=1+i; k<= i+ n[3]; k++){
		       i2 = sort2[k];
		       resid[j][i2] += xmat[j][i2] - xmean[j];
		   }
	       }
	       n[6] = n[2];
	   } else { /* the harder case, Efron approx */
	       meanwt = n[4]/n[3];  /* n[3] = number of deaths */
	       for (j=0; j<nvar; j++) xmean[j] =0;
	       for (dd=0; dd<n[3]; dd++) {
		   downwt = dd/n[3];
		   temp = n[2] - downwt* n[5];  /* working denom */
		   n[6] += 1/temp;
		   hazard = meanwt/temp;
		   cumhaz += hazard;
		   for (j=0; j<nvar; j++) {
		       tmean = (xsum1[j] - downwt*xsum2[j])/ temp;
		       xmean[j] += tmean/n[3];  /* save the average mean */
		       xhaz[j] += tmean*hazard;
		       for (k=1+i ; k<= i+ n[3]; k++) {
			   i2 = sort2[k];
			   temp2 = xmat[j][i2] - tmean;
			   resid[j][i2] += temp2/n[3];
			   resid[j][i2] += temp2 * risk[i2] *hazard *downwt;
		       }
		   }
	       }
	       n[6] = n[3]/n[6];  /* harmonic mean of the denominators */
	   }
	   /* copy per event time results to R structure */
	   rtime[itime] =  dtime;	
	   rstrat[itime] = istrat;
	   for (j=0; j<nvar; j++)
	       xbar[j][itime] = xmean[j];
	   for (j=0; j<7; j++) rn[j][itime] = n[j];
	   itime--;
  
       } /* done processing the deaths */
    } 

    /* final work for the last stratum */
    for (k= stratastart; k>=0; k--) {
	i2= sort2[k];
	for (j=0; j < nvar; j++ )
	    resid[j][i2] += risk[i2]* (xhaz[j] - xmat[j][i2]* cumhaz);
    }
 
    UNPROTECT(1);
    return(rlist);
}
    
