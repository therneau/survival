/*
** create the lagged state variable, used by the multicheck function
*/
#include "survS.h"

SEXP multicheck(SEXP y2, SEXP id2, SEXP istate2, SEXP sort2) {
    int i, ii, ilag;
    int n;
    
    double *time1, *time2, *status;
    int *istate, *id;
    int *sort;

    /* returned object*/
    SEXP pstate2;
    int *pstate;
    
    n = LENGTH(id2);
    time1 = REAL(y2);
    time2 = time1 + n;
    status = time2 + n;
    id = INTEGER(id2);
    istate = INTEGER(istate2);
    sort = INTEGER(sort2);

    PROTECT(pstate2 = allocVector(INTSXP, n));
    pstate = INTEGER(pstate2);

    ilag = sort[0];
    pstate[ilag] = istate[ilag];
	    
    for (i=1; i<n; i++) {  /* ii= next obs, ilag=prior obs */
	ii = sort[i];
	if (id[ii] == id[ilag]) {
	    if (status[ilag]>0) pstate[ii] = status[ilag];
	    else pstate[ii] = pstate[ilag];
	}	
	else pstate[ii] = istate[ii];
	ilag = ii;
    }	
    return(pstate2);
}	
    
    
