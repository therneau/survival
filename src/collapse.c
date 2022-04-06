/*
** For survival curves or their residuals, shorten a data set by
**  combining adjacent intervals such as (0, 10), (10, 15), (15,20) into
**  a single one.
** To collapse: they must have the same id, the same initial state, same wt, 
**  all but the last are not an event, and the times align.
**
** y = Surv matrix
** x = stratum
** istate = initial state
** id   = subject id
** order = order vector for sorting
*/
#include <math.h>
#include "survS.h"
#include "survproto.h"

SEXP collapse(SEXP y2,  SEXP x2, SEXP istate2, SEXP id2, SEXP wt2, 
	      SEXP order2) {
    int i, j, k, k1, k2, n;

    double *time1, *time2, *status, *wt;
    int *istate, *id, *order, *x;
    int *i1, *i2;   /* start and stop pointers */
    SEXP outmat;

    n = LENGTH(istate2);
    time1 = REAL(y2);
    time2 = time1 + n;
    status = time2 + n;
    x = INTEGER(x2);
    istate = INTEGER(istate2);
    id = INTEGER(id2);
    wt = REAL(wt2);
    order = INTEGER(order2);

    i1= (int *) R_alloc(2*n, sizeof(int));
    i2 = i1 + n;
	
    j=0; i=0;
    while (i<n) {
	k1 = order[i];
	i1[i] = k1; /* first obs for this id set */
	for (k=i+1; k<n; k++) {  /* look ahead */
	    k2 = order[k];
	    if (status[k1] !=0 || (id[k1] != id[k2]) || (x[k1] != x[k2]) ||
		(time1[k1] != time2[k2]) || (istate[k1] != istate[k2]) ||
		(wt[k1] != wt[k2])) break;
	    i++;
	    k1= k2;
	}

	i2[j] = k1; /* last obs of this set */
	i++;
	j++;
    }
	
    /* j is now the number of chains that were found */
    outmat = allocMatrix(INTSXP, j, 2);
    id = INTEGER(outmat);   /* reused id as a temp */
    for (i=0; i<j; i++) {
	id[i] = i1[i] +1;    /* +1 for R indexing */
	id[i+j] = i2[i] +1;
    }	
	
    return(outmat);
}
