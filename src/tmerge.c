/*
** Fill in a new time-dependent variable
** 
**  id: subject id in the baseline datata set (integer)
**  time1: start time for each interval in the baseline data
**  nid, ntime: the id and time point for the new covariate
**  x:  value of the new covariate
**  indx:  starting point for matching in the baseline
**
**  Both data sets are in order of time within id
*/
#include "survS.h"
#include "survproto.h"

/* First routine, for cumtdc, return a cumulative sum */
SEXP tmerge(SEXP id2,  SEXP time1x, SEXP newx2,
            SEXP nid2, SEXP ntime2, SEXP x2,  SEXP indx2) {
    int i, k;
    int n1, n2, oldid;
    int hasone =0;

    int *id, *nid;
    double *time1, 
	   *ntime, *x;
    double csum =0;
    double *newx;
    SEXP newx3;

    n1 = LENGTH(id2);  /* baseline data set */
    n2 = LENGTH(nid2); /* new data set */

    id	= INTEGER(id2);
    nid = INTEGER(nid2);
    time1 = REAL(time1x);
    ntime = REAL(ntime2);
    x     = REAL(x2);

    PROTECT(newx3 = duplicate(newx2));
    newx = REAL(newx3);

    /*
    ** i= index of baseline subject, k= index of addition row
    **  oldid = prior id, id's for baseline are integers starting with 1
    ** hasone: 0 if nothing has yet been accumlated for this subject, 1
    **  it it has
    */
    oldid = -1;  /* nobody */
    k =0;
    for (i=0; i<n1; i++) {
	if (id[i] != oldid) {
	    csum = 0;  
	    oldid = id[i];
	    hasone =0;
	}
	for (; k<n2 && nid[k] < id[i]; k++);  
	for (; k<n2 && (nid[k] == id[i]) && ntime[k] <= time1[i]; k++) {
	    csum += x[k];
	    hasone = 1;
	}
	
	if (hasone ==1) {
	    if (ISNA(newx[i])) newx[i] = csum;  /* an NA is replaced */
	    else  newx[i] = newx[i] + csum;     /* otherwise incremented */
	}
    }
    
    UNPROTECT(1);
    return(newx3);
    }

/*
** Part 2 of the code, used for tdc
**  for each row of the master data (id, time1), return the row of
**  the new data (nid, ntime2) that will provide the new data.
** Based on a last-value-carried forward rule, if the master for an id
**  had time intervals of (0,5), (5,10), (15,20) and the new data had time values
**  of  -1, 5, 11, 12, the return index would be 1, 2, and 4.  Covariates take
**  effect at the start of an interval.  Notice that obs 3 is never used.
*/
	    
SEXP tmerge2(SEXP id2,  SEXP time1x, SEXP nid2, SEXP ntime2) {
    int i, k;
    int n1, n2;

    int *id, *nid;
    double *time1, 
	   *ntime;
    SEXP index2;
    int  *index;

    n1 = LENGTH(id2);  /* baseline data set */
    n2 = LENGTH(nid2); /* new data set */

    id	= INTEGER(id2);
    nid = INTEGER(nid2);
    time1 = REAL(time1x);
    ntime = REAL(ntime2);

    PROTECT(index2 = allocVector(INTSXP, n1));
    index = INTEGER(index2);

    /*
    ** Every subject in the new data (nid, ntime) will be found in the baseline
    ** data set (id, time1) -- if not the parent routine has already tossed
    ** them, but not every obs in id will have a representative in the new.
    ** For those we return 0, otherwise the max k such that nid[k]== id[i]
    ** and ntime[k] <= time1[i]
    **
    ** For each element in (id, time1):
    **   set index to 0
    **   walk forward in data set 2 until the newid is >= current
    **   while (id matches and newtime <= oldtime), set pointer
    **        to this row
    */
    k=0;  /* index for newid */
    for (i=0; i< n1; i++) {
	index[i] =0;  /* default, assume we won't find a match */
	for (; k< n2 && (nid[k] < id[i]); k++);
	for (; k< n2 && (nid[k] == id[i]) && (ntime[k] <= time1[i]); k++) {
	    index[i] = k+1;
	}
	k--;  /* the next obs might need the same k */
    }

    UNPROTECT(1);
    return(index2);
    }


/*
** This routine is used by surv2data to accomplish last value carried
**  forward.  It started out as a modification of tmerge2, but then
**  simplified.
** The input data has id, time, and missing yes/no.  The return value
**  is a vector of integers pointing to the replacement element, or
**  0 if there is no good replacement, e.g., the first time point for a
**  subject is missing.  
** Id is an integer.
*/	    
 SEXP tmerge3(SEXP id2, SEXP miss2) {
    int i, k;
    int n;
    int oldid;

    int *id, *miss;
    SEXP index2;
    int  *index;

    n = LENGTH(id2);  /* baseline data set */
    id	= INTEGER(id2);
    miss = INTEGER(miss2);  /* actually logical, but they pass as integer*/

    PROTECT(index2 = allocVector(INTSXP, n));
    index = INTEGER(index2);

    /*
    ** The input is sorted by time within id.
    **   Simply keep track of the last non-missing row we see, resetting that
    **   constant to 0 each time a new identifier arises.
    */
    oldid = -1;   /* not anybody */
    k =0;         /* the row of interest */
    for (i=0; i<n; i++) {
	if (id[i] != oldid) {
	    k = 0;
	    oldid = id[i];
	}
	if (miss[i] ==1) index[i] =k;
	else {
	    index[i] =i +1;
	    k = i +1;      /* parent routine indices start at 1 */
	}	
    }

    UNPROTECT(1);
    return(index2);
    }

     
