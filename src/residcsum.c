/*
** A cumulative sums for each col of a matrix, resetting to 0 at each strata
**  used by the residuals.survfit routines.  If there were no strata, then
**  apply(y, 2, cumsum) does the same thing.
*/
#include "survS.h"
#include "survproto.h"
SEXP residcsum(SEXP y2, SEXP strata2)  {
    int i, j;
    int n, nc;
    double *y, temp;
    int *strata, cstrat;   /* cstrat = current strata */
    SEXP csum;
    
    PROTECT(csum= duplicate(y2));
    
    n = nrows(y2);	
    nc = ncols(y2);
    y = REAL(csum);
    strata = INTEGER(strata2);

    for (j=0; j<nc; j++) {
	for (i=0; i<n; i++) {
	    if (i==0 || strata[i] != cstrat) {
		temp =0;
		cstrat = strata[i];
	    }
	    
	    temp += *y;
	    *y++ = temp;
	}
    }
    unprotect(1);
    return(csum);
}
