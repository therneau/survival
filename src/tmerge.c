/*
** Fill in a new time-dependent variable
** 
**  id: subjec id in the baseline datata set (integer)
**  time1, time2: time intervals in the baseline data, only the second is needed
**  nid, ntime: the id and time point for the new covariate
**  x:  vaule of the new covariate
**  indx:  starting point for matching in the baseline
**
**  the baseline need not be in sorted id order, but all time points
**   for a given id are together, and in order.
*/
#include "survS.h"
SEXP tmerge(SEXP id2,  SEXP time2x, SEXP newx2,
            SEXP nid2, SEXP ntime2, SEXP x2,  SEXP indx2) {
    int i, k;
    int n1, n2;

    int *id, *nid;
    double *time2, 
	   *ntime, *x;
    int *indx;
    double *newx;
    SEXP newx3;

    n1 = LENGTH(id2);  /* baseline data set */
    n2 = LENGTH(nid2); /* new data set */

    id	= INTEGER(id2);
    nid = INTEGER(nid2);
    time2 = REAL(time2x);
    ntime = REAL(ntime2);
    x     = REAL(x2);
    indx  = INTEGER(indx2);

    PROTECT(newx3 = duplicate(newx2));
    newx = REAL(newx3);

    /*
    **	There are two indices i= baseline data set, k= additions
    **  For a given subject we might have time intervals of
    **   say (2, 5], [5,9], (12,15], (15,18]
    **   with a newtime, x pairs of (1,105), (15, 202)
    **  In this case the first 3 intervals get a value of 105 and
    **   the last a value of 202.  
    **  The indx variable says where to start for each new addition,
    **   one continues as long as time2 > newtime and id=nid.
    */
    for (k=0; k<n2; k++) {
	for (i=indx[k]-1; i<n1; i++) {
	    if (id[i] != nid[k] || time2[i] <= ntime[k]) break;
	    newx[i] = x[k];
	    }
	}
    
    UNPROTECT(1);
    return(newx3);
    }

	    
    
