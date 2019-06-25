/*
** create three variables
**   1. dupid = 1 if the id is the second or more
**   2. gap   = -1 the start time for this obs is before the prior end time
**               1 the start time for this obs is after the prior end time
**               0 times match, or the first instance of the observation
**   3. cstate = current state for this obs, obtained by chaining forward
**        initial state at entry,  first transtion, second, etc.
*/
#include "survS.h"

SEXP multicheck(SEXP y2, SEXP status2, SEXP id2, SEXP istate2, SEXP sort2) {
    int i, ii, ilag;
    int n;
    
    double *time1, *time2;
    int *status;
    int *istate, *id;
    int *sort;

    /* returned object*/
    static const char *outnames[]={"dupid", "gap", "cstate", ""};
    int *dupid, *gap, *cstate;
    SEXP rlist;   /* return list*/
    
    n = LENGTH(id2);
    time1 =  REAL(y2);
    time2 =  time1 + n;
    status = INTEGER(status2);
    id = INTEGER(id2);
    istate = INTEGER(istate2);
    sort = INTEGER(sort2);

    PROTECT(rlist = mkNamed(VECSXP, outnames));
    dupid = INTEGER(SET_VECTOR_ELT(rlist, 0, allocVector(INTSXP, n)));
    gap   = INTEGER(SET_VECTOR_ELT(rlist, 1, allocVector(INTSXP, n)));
    cstate= INTEGER(SET_VECTOR_ELT(rlist, 2, allocVector(INTSXP, n)));

    ilag = -1;
    cstate[0] = istate[0];
	    
    for (i=0; i<n; i++) {  /* ii= next obs, ilag=prior obs */
	ii = sort[i];
	if (id[ii] == id[ilag]) {
	    dupid[ii] = 1;
	    if (time1[ii] == time2[ilag]) gap[ii] =0;
	    else if (time1[ii] > time2[ilag]) gap[ii] =1;
	    else gap[ii] = -1;

	    if (status[ilag]>0) cstate[ii] = status[ilag];
	    else cstate[ii] = cstate[ilag];
	}	
	else {  /* first obs for a new id*/
	    dupid[ii] =0;
	    gap[ii] =0;
	    cstate[ii] = istate[ii];
	}	
	ilag = ii;
    }
    UNPROTECT(1);
    return(rlist);
}	
    
    
