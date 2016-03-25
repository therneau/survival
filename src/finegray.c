/*
** Do the indexing for a Fine-Gray model
**  this routine has a lot of overlap with survsplit
** input: (tstart, tstop] the time interval for each obs
**        cut: the vector of new cutpoints
**        prob: survival estimates for any new intervals
**        extend: should this observation be extended?
**        ctime, cprob: the underlying survival curve
**
** output: row = which row of the original data, for each row of the new
**         (start, end] = the new time intervals
**         wt           = probability weight for the interval.  This is 1.0
**             for all but the extended intervals.
*/
#include "survS.h"
#include "survproto.h"
#include <stdio.h>

SEXP finegray(SEXP tstart2, SEXP tstop2,   SEXP cut2,   SEXP prob2, 
	      SEXP extend2, SEXP maxtime2, SEXP ctime2, SEXP cprob2) {
    int i,j, k, extra;
    int n;   /* number of observations */
    int ncut; /* number of cuts */
    int n2;   /* number of new obs */
    int nprob;
    double tempwt;

    double *tstart, *tstop, *cut, *prob;
    double *ctime,  *cprob;
    int *extend;
    double maxtime;

    /* returned objects */
    SEXP row2, interval2, start2, end2, wt2, rlist;
    double *start, *end, *wt;
    int *row, *interval;
    static const char *outnames[]={"row", "interval", "start", "end", "wt", ""};
    n = LENGTH(tstart2);
    ncut = LENGTH(cut2);
    tstart = REAL(tstart2);
    tstop  = REAL(tstop2);
    cut   = REAL(cut2);
    prob  = REAL(prob2);
    extend= INTEGER(extend2);
    maxtime= asReal(maxtime2);
    nprob  = LENGTH(ctime2);
    ctime  = REAL(ctime2);
    cprob  = REAL(cprob2);

    /* 
    ** how many obs will I need? NA inputs are left alone.
    **  Extend observations create an extra interval for each cut beyond
    **   their endpoint.
    */
    extra =0;
    for (i=0; i<n; i++) {
	if (!ISNAN(tstart[i]) && !ISNAN(tstop[i]) ) {
	    for (j=0; j<ncut; j++) {
		if (extend[i] && cut[j] > tstop[i]) extra++;
	    }
	}	
    }

    /* allocate output */
    n2 = n + extra;
    PROTECT(rlist = mkNamed(VECSXP, outnames));
    row2 = SET_VECTOR_ELT(rlist, 0, allocVector(INTSXP, n2));
    row = INTEGER(row2);

    interval2 = SET_VECTOR_ELT(rlist, 1, allocVector(INTSXP, n2));
    interval = INTEGER(interval2);
	
    start2 = SET_VECTOR_ELT(rlist, 2, allocVector(REALSXP, n2));
    start = REAL(start2);
    end2   = SET_VECTOR_ELT(rlist, 3, allocVector(REALSXP, n2));
    end = REAL(end2);
    wt2 =    SET_VECTOR_ELT(rlist, 4, allocVector(REALSXP, n2));
    wt = REAL(wt2);
    
    /* do the work */
    k =0; 
    for (i=0; i<n; i++) {
	/* put out the original interval */
	start[k] = tstart[i];
	end[k]   = tstop[i];
	row[k]  = i+1;    /* 1 based subscripts for R */
	interval[k] = 1;
	wt[k] =1;
	k++;
    
	if (extend[i] && tstop[i] < maxtime) {
	    /* Add more itervals at the end */
	    /* First skip over cuts that we don't need */
	    tempwt = 1;
	    for (j=0; (j<nprob) && (ctime[j] <= tstop[i]); j++) {
		tempwt = cprob[j];
		printf("i= %d, j=%d, tempwt=%f\n", i, j, tempwt);
	    }		
	    for (j=0; (j<ncut) && (cut[j] <= tstop[i]); j++);
		
	    start[k] = tstop[i];
	    row[k] = i+1;
	    interval[k] = j;
	    
	    for (; j<ncut; j++) {
		end[k] = cut[j];
		wt[k] = prob[j]/tempwt;
		k++;
		start[k] = cut[j];
		row[k] =  i+1;
		interval[k] = j+1;
		}
	    wt[k] = prob[j]/tempwt;
	    end[k] = maxtime;
	    k++;
	}
    }
    UNPROTECT(1);
    return(rlist);
    }
