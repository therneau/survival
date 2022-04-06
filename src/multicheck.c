/*
** create three variables
**   1. dupid = 1 first obs for this id, 2=last, 3=both, 0=neither
**   2. gap   = -1 the start time for this obs is before the prior end time
**               1 the start time for this obs is after the prior end time
**               0 times match, or the first instance of the observation
**   3. cstate = current state for this obs, obtained by chaining forward
**        initial state at entry,  first transtion, second, etc. but ignoring
**	  censors as a 'change'.
*/
#include "survS.h"
#include "survproto.h"

SEXP multicheck(SEXP time12, SEXP time22, SEXP status2, SEXP id2, 
		SEXP istate2, SEXP sort2) {
    int i, ii, oldid, oldii;
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
    time1 =  REAL(time12);
    time2 =  REAL(time22);
    status = INTEGER(status2);
    id = INTEGER(id2);
    istate = INTEGER(istate2);
    sort = INTEGER(sort2);

    PROTECT(rlist = mkNamed(VECSXP, outnames));
    dupid = INTEGER(SET_VECTOR_ELT(rlist, 0, allocVector(INTSXP, n)));
    gap   = INTEGER(SET_VECTOR_ELT(rlist, 1, allocVector(INTSXP, n)));
    cstate= INTEGER(SET_VECTOR_ELT(rlist, 2, allocVector(INTSXP, n)));

    oldid = -1;  /* the input id values are all > 0 */
    oldii = 0 ;  /* stop a complaint from -Wall option of gcc */
    cstate[0] = istate[0];
	    
    for (i=0; i<n; i++) {  /* ii= next obs, ilag=prior obs */
	ii = sort[i];
	if (id[ii] == oldid) {
	    dupid[ii] = 0;
	    if (time1[ii] == time2[oldii]) gap[ii] =0;
	    else if (time1[ii] > time2[oldii]) gap[ii] =1;
	    else gap[ii] = -1;

	    if (status[oldii]>0) cstate[ii] = status[oldii];
	    else cstate[ii] = cstate[oldii];
	}	
	else {  /* first obs for a new id*/
	    oldid = id[ii];
	    dupid[ii] =0;
	    gap[ii] =0;
	    cstate[ii] = istate[ii];
	    if (i>0) dupid[oldii] += 2;  /* prior obs was last for that id*/
	}	
	oldii = ii;
    }
    dupid[oldii] += 2;
    UNPROTECT(1);
    return(rlist);
}	
    
    
