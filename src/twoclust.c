/*
** a check called by survfitAJ: does any id appear in multiple clusters?
** Input:
**  id: integer vector of id values
**  cluster: integer vector if cluster values
**  idord: ordering vector for the id value
**  Return: 0 = no, 1= yes
*/
#include "survS.h"

SEXP twoclust(SEXP id2,  SEXP cluster2,  SEXP idord2) {
    int i, n, cid;
    int iclust;  /* cluster to which the current id belongs */
    int *id, *cluster, *idord;
    SEXP rval2;
    int *rval;


    /* create the output object */
    rval2 = PROTECT(allocVector(INTSXP, 1));
    rval = INTEGER(rval2);

    n = length(id2);
    id = INTEGER(id2);
    cluster = INTEGER(cluster2);
    idord  = INTEGER(idord2);

    for (i=0; i<n; ) { 
	cid = id[idord[i]];  /* current id */
	iclust = cluster[idord[i]];
	for (; i<n && id[idord[i]] == cid; i++) {
	    if (cluster[idord[i]] != iclust) {
		*rval =1;
		UNPROTECT(1);
		return(rval2);
	    }
	}
    }
	    
    *rval =0; 
    UNPROTECT(1);
    return(rval2);
}
