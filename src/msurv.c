/* Automatically generated from all.nw using noweb */
/*  -*- c -*-  */
#include "Rinternals.h"
SEXP msurv(SEXP nrisk2, SEXP wrisk2, SEXP nevent2, SEXP wevent2,
           SEXP itime2, SEXP status2,
           SEXP prior2, SEXP etype2,  SEXP wt2,
           SEXP dsort2, SEXP isort2) {
    int i,j,k;
    double *wrisk, *wevent, *wt;
    int    *nrisk, *status, *prior, *etype, *nevent;
    int    *dtime,  *etime;  /* death time and optional entry time */
    int    *dsort, *isort;
    int    eflag, i2, k2;
    SEXP   retlist;  
    int    time;  /* current time of interest */
    int    ntime, nstate, n;
    static const char *outnames[]= {"nrisk", "wrisk", "nevent", "wevent", ""};
    
    nrisk = INTEGER(nrisk2);
    wrisk = REAL(wrisk2);
    wevent= REAL(wevent2);
    nevent= INTEGER(nevent2);
    dtime = INTEGER(itime2);
    status= INTEGER(status2);
    prior = INTEGER(prior2);
    etype = INTEGER(etype2);
    wt    = REAL(wt2);
    dsort = INTEGER(dsort2);
    nstate= nrows(nrisk2);
    n     = length(dsort2);
    
    /* 
    ** Add up the risk set and the deaths 
    **   Walk backwards through the observations and through time
    */
    if (ncols(itime2)==2) {
        etime = dtime;
        dtime += n;  /* point at the death time */
        isort = INTEGER(isort2);
        eflag=1;
        i2 = n-1;
        k2 = isort[i2];
        }
    else eflag=0;
    
    for (i=n-1; i>=0; ) {
        k = dsort[i];
        time = dtime[k];  /* current time of interest (there may be ties) */

        while (eflag==1 && i2>=0 && etime[k2] >=time) {
            /* remove those who start later than "time" from the risk set */
            wrisk[prior[k2]] -= wt[k2];
            nrisk[prior[k2]] --;
            i2--;
            k2 = isort[i2];
        }
 
        
        if (i<(n-1)) {
            /* new death time */
            for (j=0; j<nstate; j++) nrisk[j+nstate] = nrisk[j];
            nrisk += nstate;
            wrisk += nstate;
            nevent += nstate*nstate;
            wevent += nstate*nstate;
            }
            
        while(i >=0 && dtime[k] == time) {
            if (status[k] ==1) {
                nevent[prior[k]+ nstate*etype[k]]++;
                wevent[prior[k] +nstate*etype[k]] += wt[k];
            }
            wrisk[prior[k]] += wt[k];
            nrisk[prior[k]] ++;
            i--;
            k = dsort[i];
            }
        }
    
    /*
    ** Create output structure
    */
    PROTECT(retlist = mkNamed(VECSXP, outnames));
    SET_VECTOR_ELT(retlist, 0, nrisk2);
    SET_VECTOR_ELT(retlist, 1, wrisk2);
    SET_VECTOR_ELT(retlist, 2, nevent2);
    SET_VECTOR_ELT(retlist, 3, wevent2);
    UNPROTECT(1);
    return(retlist);
}        
