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
** The input data is assumed sorted by strata, decreasing y, and censors before
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
** This is the (time1, time2, status) version.
** y:     the survival times, in data set order
** sort1: ordering vector: strata, decreasing time1
** sort2: ordering vector, strata, decreasing time2, status
** strat: will be 1 if sort2[i] is the first obs in a new stratum
**   
** The hard part of this routine is that we want to keep a list of all who
** are at risk, other routines just need sums.
** An easy method is to keep an atrisk[] vector, in the same order as y, whose
** elements are set to 1 when someone enters the risk set (they appear in sort2)
** and returns to 0 when they leave (appear in sort1).  For each unique death
** time scan over all n elements of the vector and write them out.  At the
** start of each strata make sure to set them all to zero.  If there are lots
** of strata then the average size of a risk set will be much less than n, and
** this will be very inefficient.
**  
** However, for a large data set the tt() intermediate data is huge, whatever we
** do.  People should not be using tt() in that case.  That won't stop users from
** doing so, of course.
**
** This routine keeps two ancillary vectors to modestly speed this up.  The
** who vector has the row numbers of those at risk, in elements 0:(nrisk-1)
** The second vector atrisk is a reverse index: if who[i] = k, then atrisk[k] =i.
** When someone is added to the set of subjects at risk, add them to the end
** of who[] and update atrisk[].
** When someone is removed, first look them up in atrisk[], say they are element 
**  k in who
**   - set who[k] = who[nrisk--]  (the order of elements in who is immaterial)
**   - update the element of atrisk that just moved
**   - we don't have to set the old element of atrisk to 0 since it will never
**     be looked at again.
*/

SEXP coxcount2(SEXP y2, SEXP isort1, SEXP isort2, SEXP strat2) {
    int ntime, nrow;
    int i, j, k, n;
    int nrisk=0, *atrisk, *who;
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
    ** Be awake to an odd special case: the current strata finishes with
    **  an event at time 25, say, and the next strata starts with an event
    **  at time 25.  They need to be treated as separate times.  For the
    **  removal step below "j<i" suffices -- the risk set can't have 0 subjects.
    **  For the addition, don't intrude into a new stratum
    */
    ntime=0; nrow=0;
    j =0;  /* walks along the sort1 vector (start times) */
    for (i=0; i<n;) {
	iptr = sort2[i];
        if (strata[i]== 1) {
	    nrisk=0; 
	    j = i;    /* reset */
	}
	
        if (status[iptr] ==1) {
            ntime++;
            dtime = time2[iptr];
	    /* remove those who start after this death */
            for (; j<i && time1[sort1[j]] >= dtime; j++) 
                         nrisk--;
	    /* add this death, then any ties */
	    nrisk++; i++;  

            for(; i< n && strata[sort2[i]]==0 && time2[sort2[i]] == dtime; i++) {
		nrisk++;
	    }
            nrow += nrisk;
        } else {
	    nrisk++;
	    i++;
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
    atrisk = (int *)R_alloc(2*n, sizeof(int)); /* marks who is at risk */
    who    = atrisk + n;
    
    /*
    ** Pass 2, fill them in
    */
    ntime=0; nrisk=0;
    j=0;  /* pointer to time1 */;
    for (i=0; i<n; ) {
        iptr = sort2[i];
        if (strata[i] ==1) { 
	    nrisk =0;
	    j =i;    /* start anew on removals too */
	}

	/* add this person to the risk set */
	if (status[iptr] ==0) {
	    atrisk[iptr] = nrisk;
	    who[nrisk++] = iptr;
 	    i++;
	}
	else {
            dtime = time2[iptr];
	    /* unmark those who are no longer at risk */
            for (; j <i && time1[sort1[j]] >=dtime; j++) {
		jptr = sort1[j];    /* who should be removed */
		k = atrisk[jptr];   /* is stored here in who[] */
		who[k] = who[--nrisk]; /* someone else takes their place */
		atrisk[who[k]] = k;    /* update pointer for the usurper */
	    }

	    /* write out the index & status for the controls */
            for (k=0; k< nrisk; k++) {
		*rstatus++ =0;
		*rindex++  = who[k] +1;  /* R subscripts start at 1 */
	    }
	    
	    /* add this death */
	    *rstatus++ =1;
	    *rindex++  = iptr +1;
	    atrisk[iptr] = nrisk;
	    who[nrisk++] = iptr;
	    i++;
	    /* walk through any tied deaths */
	    for (; i<n && strata[sort2[i]]==0 && time2[sort2[i]]==dtime; i++) {
		iptr = sort2[i];
		*rstatus++ = 1; 
		*rindex++ = iptr +1 ; /* R subscripts start at 1 */
		atrisk[iptr] = nrisk;
		who[nrisk++] = iptr;
	    }
 
            rtime[ntime] = dtime;
            rn[ntime] = nrisk;
            ntime++;
        }
    }    

    UNPROTECT(1);
    return(rlist);
}
