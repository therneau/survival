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
	      SEXP extend2, SEXP ctime2, SEXP cprob2) {
    int i,j, k, extra;
    int n;   /* number of observations */
    int ncut; /* number of cuts */
    int n2;   /* number of new obs */
    int nprob;
    double tempwt;

    double *tstart, *tstop, *cut;
    double *ctime,  *cprob, *prob;
    int *extend;

    /* returned objects */
    SEXP row2, start2, end2, wt2, rlist;
    double *start, *end, *wt;
    int *row;
    static const char *outnames[]={"row", "start", "end", "wt", ""};

    n = LENGTH(tstart2);
    ncut = LENGTH(cut2);
    tstart = REAL(tstart2);
    tstop  = REAL(tstop2);
    cut   = REAL(cut2);
    prob  = REAL(prob2);
    extend= INTEGER(extend2);
    nprob  = LENGTH(ctime2);
    ctime  = REAL(ctime2);
    cprob  = REAL(cprob2);

    /* 
    ** how many obs will I need? NA inputs are left alone.
    **  Extend observations have weight 1 up to the next cutpoint after their
    **  max, and an extra for any cutpoints after that.
    */
    extra =0;
    for (i=0; i<n; i++) {
	if (!ISNAN(tstart[i]) && !ISNAN(tstop[i]) && extend[i] ) {
	    for (j=0; (j<ncut) && (cut[j] <= tstop[i]); j++);
	    if (j < ncut) extra += (ncut -j);
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
    
    /* do the work */
    k =0; 
    for (i=0; i<n; i++) {
	/* put out the original interval */
	start[k] = tstart[i];
	end[k]   = tstop[i];
	row[k]  = i+1;    /* 1 based subscripts for R */
	wt[k] =1;
    
	if (extend[i] && cut[ncut-1] > tstop[i]) {
	    /* 
	    ** Add more itervals at the end 
	    **  tempwt = weight for the interval in which they ended
	    **  If they are removed at exactly a cutpoint, the assumption
	    **  is that censoring was not possible on that given day.
	    */
	    tempwt = 1;
	    for (j=0; (j<nprob) && (ctime[j] < tstop[i]); j++) {
		tempwt = cprob[j];
	    }		
	    printf("i= %d, stop=%f, tempwt=%f\n", i, tstop[i], tempwt);

	    for (j=0; (j<ncut) && (cut[j] <= tstop[i]); j++);
	    printf(" j=%d, cut=%f, stop=%f\n", j, cut[j], tstop[i]);
	    if (prob[j] == tempwt) {
		/* tstop < cut, but the censoring distribution has no
		   jumps before that cut
		*/
		end[k] = cut[j];
		j++;
		}
	    if (j < ncut) {
		k++;
		row[k] = i+1;
		start[k] = tstop[i];
		end[k]  = cut[j];
		wt[k]   = prob[j]/tempwt;
		for (j=j+1; j<ncut; j++) {
		    k++;
		    row[k] = i+1;
		    start[k] = cut[j-1];
		    end[k] = cut[j];
		    wt[k] = prob[j]/tempwt;
		}
	    }
	}
	k++;
    }
    UNPROTECT(1);
    return(rlist);
}
