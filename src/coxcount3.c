/*
** Counts for the survfit.coxphms routine
**
** y       : 3 column of time1, time2, 0/1 
** weight  : case weight
** rwt     : risk weight
** sort1   : order vector for time1
** sort2   : order vector for time2
** strata  : strata 1, 2, 3 etc
** transition: which transition
** ntran   : number of hazard functions (some may be shared)
** utime   : the times at which to report the counts
*/

#include "survS.h"
#include "survproto.h"
#include <math.h>

SEXP coxcount3(SEXP y2,      SEXP weight2,    SEXP rwt2,
	       SEXP sort12,  SEXP sort22,     SEXP strata2,
	       SEXP transition2, SEXP nhaz2,  SEXP utime2) {
    
    /* pointers to the R input variables */
    int *sort1, *sort2;  /*sort indexs for time1 and time2 */
    double *time1, 
	   *time2,
	   *status;     /* (time1, time2, state) survival time */
    double *wt, *rwt;   
    int *strata, *transition;
    double *utime;
    int nhaz;
    int ntime;
    int nused;

    /* returned and working objects */
    SEXP rlist;         /* the returned list and variable names of same */  
    const char *rnames[]= {"nrisk", "nevent", "cumhaz", ""};
    SEXP setemp;
    double **nrisk, **nevent, **cumhaz;
    
    /* working variables */
    int i, j, k, dindex, first;
    int i1, i2;   /* index the start/stop variables */
    int istrat;   /* current strata */
    double dtime; /* current event time */
    int debug;
    
    nused = LENGTH(weight2);
    time1 = REAL(y2);
    time2 = time1 + nused;
    status= time2 + nused;
    wt    = REAL(weight2);
    rwt   = REAL(rwt2);
    sort1 = INTEGER(sort12);
    sort2 = INTEGER(sort22);
    strata= INTEGER(strata2);
    transition = INTEGER(transition2);
    nhaz =  asInteger(nhaz2);
    utime = REAL(utime2);
    ntime = LENGTH(utime2);
       
    PROTECT(rlist=mkNamed(VECSXP, rnames));
    setemp = SET_VECTOR_ELT(rlist, 0, allocMatrix(REALSXP, ntime, nhaz));
    nrisk  =  dmatrix(REAL(setemp), ntime, nhaz);
    setemp = SET_VECTOR_ELT(rlist, 1, allocMatrix(REALSXP, ntime, nhaz));
    nevent  = dmatrix(REAL(setemp), ntime, nhaz);
    setemp = SET_VECTOR_ELT(rlist, 2,  allocMatrix(REALSXP, ntime, nhaz));
    cumhaz  = dmatrix(REAL(setemp), ntime, nhaz);
 
    /* zero the return values */
    for (i=0; i<ntime; i++) {
	for (j=0; j< nhaz; j++) {
	    nrisk[j][i] =0;
	    nevent[j][i] =0;
	    cumhaz[j][i] =0;
	}
    }
    
    /*
    ** Because an observation can have different rwt values for different
    **  transitions, the input data has a row for every subject and every
    **  transition for which they are at risk. 
    ** The output has a row for each unique event time (per strata)
    ** Walk backwards through the data 
    */
    i1 = nused -1;
    i2 = nused -1; /* add and remove from the risk set */
    dindex = ntime-1;
    istrat= strata[nused -1];
    debug=0;
    for (i= nused-1; i>=0; i--) {
	debug++;
	if (debug< 0)printf("i= %d, istrat=%d, sort=%d,strata=%d\n", i, istrat,
	       sort2[i], strata[sort2[i]]); 
	if (i == (nused -1) || istrat!= strata[sort2[i]]) first=1; 
	istrat = strata[i]; /* current strata */
	dtime = utime[dindex];
	if (debug <0) printf("debug=%d, dindex=%d, dtime=%f\n", debug, dindex, dtime); 

	/* remove any who have left the risk set */
	for (; i1>=0 && time1[sort1[i1]] >= dtime && strata[i1]== istrat; i1--){
	    j = sort1[i1];
	    k = transition[j];
	    nrisk[k][dindex] = nrisk[k][dindex] - wt[j]*rwt[j];
	}
	    
	/* add into the risk set */
	if (debug<0) printf("dtime=%f, dindex=%i, i2=%d, sort=%d, time2=%f\n",
			    dtime, dindex, i2, sort2[i2], time2[sort2[i2]]);
	for (; i2>=0 && time2[sort2[i2]]>= dtime; i2--) {
	    j = sort2[i2];
	    k = transition[j] -1L;
	if (debug<0) printf("i2=%d, sort=%d, time2=%f k=%d\n",
			    i2,  j, time2[j], k);
	    nrisk[k][dindex]  = nrisk[k][dindex]  + wt[j]*rwt[j];
	    nevent[k][dindex] = nevent[k][dindex] + wt[j]*status[j];
	}
	/* update cumhaz */
	if (first==1) {
	    first =0;
	    for (j=0; j< nhaz; j++) 
		cumhaz[j][dindex]= nevent[j][dindex]/nrisk[j][dindex];
	} else {
	   for (j=0; j< nhaz; j++) 
	    	cumhaz[j][dindex]= nevent[j][dindex]/nrisk[j][dindex] +
		    cumhaz[j][dindex+1];
	}
	dindex--;
	i = i2;
    }
    
    UNPROTECT(1);
    return(rlist);
}
