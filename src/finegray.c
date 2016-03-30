/*
** Do the indexing for a Fine-Gray model
** input: (tstart, tstop] the time interval for each obs
**        ctime, cprob: the underlying survival curve
**        the interval that ends at ctime has probability cprob
**        extend: should this observation be extended?
**        keep: output the interval from ctime[i] to ctime[i+1]?
**
** output: row = which row of the original data, for each row of the new
**         (start, end] = the new time intervals
**          wt          = probability weight for the interval.  This is 1.0
**             for all but the extended intervals.
*/
#include "survS.h"
#include "survproto.h"
#include <stdio.h>

SEXP finegray(SEXP tstart2, SEXP tstop2,   SEXP ctime2,   SEXP cprob2, 
	      SEXP extend2, SEXP keep2) {
    int i,j, k, iadd, extra;
    int n;   /* number of observations */
    int ncut; /* number of censoring intervals */
    int n2;   /* number of new obs */
    double tempwt;

    double *tstart, *tstop;
    double *ctime,  *cprob;
    int *extend, *keep;

    /* returned objects */
    SEXP row2, start2, end2, wt2, add2, rlist;
    double *start, *end, *wt;
    int *row, *add;
    static const char *outnames[]={"row", "start", "end", "wt", "add",  ""};

    n = LENGTH(tstart2);
    ncut = LENGTH(cprob2);
    tstart = REAL(tstart2);
    tstop  = REAL(tstop2);
    extend= LOGICAL(extend2);
    keep  = LOGICAL(keep2);
    ctime = REAL(ctime2);
    cprob = REAL(cprob2);

    /* 
    ** how many obs will I need? NA inputs are left alone.
    **  Extend observations have weight 1 up to the next cutpoint after their
    **  max, and an extra for any cutpoints after that.
    */
    extra =0;
    for (i=0; i<n; i++) {
	if (!ISNAN(tstart[i]) && !ISNAN(tstop[i]) && extend[i] ) {
	    for (j=0; (j<ncut) && (ctime[j] < tstop[i]); j++);
	    /* j = the interval they lie in, 0 = before the first cut
	       point, 1 = within the first interval.  We have to add 
	       any intervals after the one the subject falls in;
	    */
	    j++;
	    for (; j<ncut; j++) extra += keep[j];
	}
    }	

    /* allocate output */
    n2 = n + extra;
    PROTECT(rlist = mkNamed(VECSXP, outnames));
    row2 = SET_VECTOR_ELT(rlist, 0, allocVector(INTSXP, n2));
    row = INTEGER(row2);

    start2 = SET_VECTOR_ELT(rlist, 1, allocVector(REALSXP, n2));
    start = REAL(start2);
    end2   = SET_VECTOR_ELT(rlist, 2, allocVector(REALSXP, n2));
    end = REAL(end2);
    wt2 =    SET_VECTOR_ELT(rlist, 3, allocVector(REALSXP, n2));
    wt = REAL(wt2);
    add2 =   SET_VECTOR_ELT(rlist, 4, allocVector(INTSXP, n2));
    add = integer(add2);
    
    /* do the work */
    k =0; 
    for (i=0; i<n; i++) {
	/* put out the original interval */
	start[k] = tstart[i];
	end[k]   = tstop[i];
	row[k]  = i+1;    /* 1 based subscripts for R */
	wt[k] =1;
	add[k] =0;
    
	if (!ISNAN(tstart[i]) && !ISNAN(tstop[i]) && extend[i]) {
	    /* ctime contains the time at the end of the interval */
	    for (j=0; (j<ncut) && (ctime[j] < tstop[i]); j++);
	    end[k] = ctime[j];  /* extend them to the end of the interval */
	    tempwt = cprob[j];

	    j++;
	    iadd=0;
	    for (; j < ncut; j++) {
		if (keep[j]) {
		    /* add more */
		    k++;
		    iadd++;
		    row[k] = i+1;
		    start[k] = ctime[j-1];
		    end[k]  = ctime[j];
		    wt[k]   = cprob[j]/tempwt;
		    add[k]  = iadd;
		}
	    }
	}	
	k++;
    }
    UNPROTECT(1);
    return(rlist);
}
