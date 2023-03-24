/*
** This routine is used by coxph when there is a tt() term in the model
**   In that case, the data set gets expanded to long form, which has a 
** separate set of rows for each unique death time.  Each event time and its
** risk set is treated as a separte stratum.  We might expect a data
** set with 1000 rows and 200 deaths to end up with 200*500 rows, if the deaths
** were spread out evenly over time.  These expansions are done separately
** within each strata.
**
** This routine returns vectors with the unique event times, the number of obs
** in the risk set at each time (ntime elements), the row indices and the status
** The latter two have the values for strata 1, then strata 2, ... (nrow
** elements).  This process is not particularly fast: it is just brute force.
**
** The inoput data is assumed sorted by strata, decreasing y, and censors before
** deaths.  The strat variable is =1 for the first obs of each stratum.
*/
 
#include "survS.h"
/*
** Count up risk sets and identify who is in each
*/
SEXP coxcount1(SEXP y2, SEXP strat2) {
    int ntime, nrow;
    int i, j, n;
    int stratastart=0;  /* start row for this strata */
    int nrisk=0;  /* number at risk (=0 to stop -Wall complaint)*/
    double *time, *status;
    int *strata;
    double dtime;
    SEXP rlist;
    const char *rnames[] = {"time", "nrisk", "index", "status", ""}; 
    double *rtime;
    int *rn, *rindex, *rstatus;
    
    n = nrows(y2);
    time = REAL(y2);
    status = time +n;
    strata = INTEGER(strat2);
    
    /* 
    ** First pass: count the total number of death times (risk sets)
    **  and the total number of rows in the new data set.
    */
    ntime=0; nrow=0;
    for (i=0; i<n; i++) {
        if (strata[i] ==1) nrisk =0;
        nrisk++;
        if (status[i] ==1) {
            ntime++;
            dtime = time[i];
            /* walk across tied times, if any */
            for (j=i+1; j<n && time[j]==dtime && status[j]==1 && strata[j]==0;
                 j++) nrisk++;
            i = j-1;
            nrow += nrisk;
        }
    }
    /*
    **  Allocate memory
    */
    PROTECT(rlist = mkNamed(VECSXP, rnames));
    rtime  = REAL(SET_VECTOR_ELT(rlist, 0, allocVector(REALSXP, ntime)));
    rn     = INTEGER(SET_VECTOR_ELT(rlist, 1, allocVector(INTSXP, ntime)));
    rindex = INTEGER(SET_VECTOR_ELT(rlist, 2, allocVector(INTSXP, nrow)));
    rstatus= INTEGER(SET_VECTOR_ELT(rlist, 3, allocVector(INTSXP, nrow)));
    
    /*
    ** Pass 2, fill them in
    */
    ntime=0; 
    for (i=0; i<n; i++) {
        if (strata[i] ==1) stratastart =i;
        if (status[i]==1) {
            dtime = time[i];
            for (j=stratastart; j<i; j++) *rstatus++=0; /*non-deaths */
            *rstatus++ =1; /* this death */
            /* tied deaths */
            for(j= i+1; j<n && status[j]==1 && time[j]==dtime  && strata[j]==0;
                j++) *rstatus++ =1;
            i = j-1;

            rtime[ntime] = dtime;
            rn[ntime] = i +1 -stratastart;
            ntime++;
            for (j=stratastart; j<=i; j++) *rindex++ = j+1;
	}
    }

    unprotect(1);
    return(rlist);
}

/*
** This is the (time1, time2, status) version
*/

SEXP coxcount2(SEXP y2, SEXP isort1, SEXP isort2, SEXP strat2) {
    int ntime, nrow;
    int i, j, istart, n;
    int nrisk=0, *atrisk;
    double *time1, *time2, *status;
    int *strata;
    double dtime;
    int iptr, jptr;

    SEXP rlist;
    const char *rnames[] = {"time", "nrisk", "index", "status", ""}; 
    double *rtime;
    int *rn, *rindex, *rstatus;
    int *sort1, *sort2;
    
    n = nrows(y2);
    time1 = REAL(y2);
    time2 =  time1+n;
    status = time2 +n;
    strata = INTEGER(strat2);
    sort1 = INTEGER(isort1);
    sort2 = INTEGER(isort2);
    
    /* 
    ** First pass: count the total number of death times (risk sets)
    **  and the total number of rows in the new data set
    */
    ntime=0; nrow=0;
    istart =0;  /* walks along the sort1 vector (start times) */
        for (i=0; i<n; i++) {
        iptr = sort2[i];
        if (strata[iptr]==1) {
	    nrisk=0;
	    istart =i;
	    }
 
	nrisk++;
        if (status[iptr] ==1) {
            ntime++;
            dtime = time2[iptr];
            for (; istart <i && time1[sort1[istart]] >= dtime; istart++) 
                         nrisk--;
            for(j= i+1; j<n; j++) {
                jptr = sort2[j];
                if (status[jptr]==1 && time2[jptr]==dtime && strata[jptr]==0)
                    nrisk++;
                else break;
                }
            i= j-1;
            nrow += nrisk;
            }
        }

    /*
    **  Allocate memory
    */
    PROTECT(rlist = mkNamed(VECSXP, rnames));
    rtime  = REAL(SET_VECTOR_ELT(rlist, 0, allocVector(REALSXP, ntime)));
    rn     = INTEGER(SET_VECTOR_ELT(rlist, 1, allocVector(INTSXP, ntime)));
    rindex = INTEGER(SET_VECTOR_ELT(rlist, 2, allocVector(INTSXP, nrow)));
    rstatus= INTEGER(SET_VECTOR_ELT(rlist, 3, allocVector(INTSXP, nrow)));
    atrisk = (int *)R_alloc(n, sizeof(int)); /* marks who is at risk */
    
    /*
    ** Pass 2, fill them in
    */
    ntime=0; nrisk=0;
    j=0;  /* pointer to time1 */;
    istart=0;
    for (i=0; i<n; ) {
        iptr = sort2[i];
        if (strata[i] ==1) {
            nrisk=0;
	    istart =i;
            for (j=0; j<n; j++) atrisk[j] =0;
            }
        nrisk++;
        if (status[iptr]==1) {
            dtime = time2[iptr];
            for (; istart<i && time1[sort1[istart]] >=dtime; istart++) {
                atrisk[sort1[istart]]=0;
                nrisk--;
                }
            for (j=1; j<nrisk; j++) *rstatus++ =0;
            for (j=0; j<n; j++) if (atrisk[j]) *rindex++ = j+1; 

            atrisk[iptr] =1;
            *rstatus++ =1; 
            *rindex++ = iptr +1; /* R subscripts start at 1 */
            for (j=i+1; j<n; j++) {
                jptr = sort2[j];
                if (time2[jptr]==dtime && status[jptr]==1 && strata[jptr]==0){
                    atrisk[jptr] =1;
                    *rstatus++ =1;
                    *rindex++ = jptr +1;
                    nrisk++;
                    }
                else break;
                }
            i = j;
            rtime[ntime] = dtime;
            rn[ntime] = nrisk;
            ntime++;
        }
        else {
            atrisk[iptr] =1;
            i++;
        }
    }    

    unprotect(1);
    return(rlist);
}
