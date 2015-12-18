/*
** Do the indexing for survSplit
** input: (tstart, tstop] the time interval for each obs
**        cut: the vector of new cutpoints
** output: row = which row of the original data, for each row of the new
**         (start, end] = the new time intervals
**         censor       = the new should be censored (new end point)
*/
#include "survS.h"
#include "survproto.h"
#include <stdio.h>

SEXP survsplit(SEXP tstart2, SEXP tstop2, SEXP cut2) {
    int i,j, k, extra;
    int n;   /* number of observations */
    int ncut; /* number of cuts */
    int n2;   /* number of new obs */

    double *tstart, *tstop, *cut;

    /* returned objects */
    SEXP row2, interval2, start2, end2, censor2, rlist;
    double *start, *end;
    int *row, *censor, *interval;
    static const char *outnames[]={"row", "interval", "start", "end", 
				   "censor", ""};
    n = LENGTH(tstart2);
    ncut = LENGTH(cut2);
    tstart = REAL(tstart2);
    tstop  = REAL(tstop2);
    cut   = REAL(cut2);

    /* 
    ** how many obs will I need?  Each cutpoint strictly within an interal
    **  generates an extra line.  NA inputs are left alone.
    */
    extra =0;
    for (i=0; i<n; i++) {
	for (j=0; j<ncut; j++)
	    if (!ISNAN(tstart[i]) && !ISNAN(tstop[i]) &&
			cut[j] > tstart[i] && cut[j] < tstop[i]) extra++;
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
    censor2= SET_VECTOR_ELT(rlist, 4, allocVector(LGLSXP,  n2));
    censor = LOGICAL(censor2);
    
    /* do the work */
    k =0; 
    for (i=0; i<n; i++) {
	if (ISNAN(tstart[i]) || ISNAN(tstop[i])) {
	    start[k] = tstart[i];
	    end[k]   = tstop[i];
	    row[k]  = i+1;    /* 1 based subscripts for R */
	    interval[k] = 1;
	    k++;
	    }
	else {
	    for (j=0; j<ncut & (cut[j]<= tstart[i]); j++);
	    start[k] = tstart[i];
	    row[k] =i+1;
	    interval[k] = j;
	    for (; j<ncut & (cut[j] < tstop[i]); j++) {
		if (cut[j] > tstart[i]) {
		    end[k] = cut[j];
		    censor[k] =1;
		    k++;
		    start[k] = cut[j];
		    row[k] =i+1;
		    interval[k] = j+1;
		}
	    }
	    end[k] = tstop[i];
	    censor[k] =0;
	    k++;
	}
    }
    UNPROTECT(1);
    return(rlist);
    }

