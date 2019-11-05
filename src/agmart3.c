/*
** The fastest possible calculation (at least the fastest I've come up with)
** for martingale residuals from a (start, stop] model.  
**
** Input
**      surv    start, stop, event matrix of survival data
**      score   the vector of subject scores, i.e., exp(beta*z)
**      weight  case weights
**      strata  integer vector of strata identifiers
**      sort     sort order for the obs: col1 = stop times, 2 = start times
**         time within strata for both
**      method  will be ==1 for the Efron method
** Output
**      martingale residual
** 
** If an observation is never at risk, e.g. (50, 72] and there are no events
**   in that window, then ensure that we never add/subtract that obs from
**   any totals.  Sometimes users "dice" a data set into a lot of small
**   intervals, and this step helps preserve numeric accuracy.  Such obs
**   have a residual of zero.
*/
#include <math.h>
#include "survS.h"
#include "survproto.h"

SEXP agmart3(SEXP nused2,  SEXP surv2,  SEXP score2, SEXP weight2, 
	     SEXP strata2, SEXP sort12, SEXP sort22, SEXP method2) {

    int i, k;
    int p1, p2, istrat;
    double deaths, denom, e_denom;
    double hazard, e_hazard, cumhaz;
    double temp, dtime =0;  /* =0 to stop a -Wall message */
    double wtsum;
    int nr, person1, person2, nused;
 
    /* pointers to the input data */
    double *tstart, *tstop, *event;
    double *weight, *score;
    int    *sort1, *sort2, *strata;

    int method;  /* integer version of input */
    int *atrisk; /* 1= ever at risk */

    /* output */
    SEXP resid2;
    double *resid;

    nused = asInteger(nused2);
    nr = nrows(surv2);
    method = asInteger(method2);
    tstart = REAL(surv2);
    tstop  = tstart +nr;
    event = tstop +nr;
    weight= REAL(weight2);
    score = REAL(score2);
    sort1 = INTEGER(sort12);
    sort2 = INTEGER(sort22);
    strata= INTEGER(strata2);

    PROTECT(resid2 = allocVector(REALSXP, nr));
    resid = REAL(resid2);
    atrisk = (int *) ALLOC(nr, sizeof(int));
    for (k=0; k<nr; k++) {
	resid[k] =0;
	atrisk[k] =0;
    }	
    /*
    **  person1, p1, sort1 refer to start times, and person2, p2, sort2
    **    to the ending times for each subject.
    **  'person2' walks through the the data from 1 to nused,
    **     sort2[0] points to the largest stop time, sort2[1] the next, ...
    **  'time' is a scratch variable holding the time of current interest
    **  'person1' walks through the start times.  It will be smaller than 
    **    'person2': if person2=27 that means that 27 subjects have stop >=time,
    **    and are thus potential members of the risk set.  If 'person1' =9,
    **    that means that 9 subjects have start >=time and thus are NOT part
    **    of the risk set.  (stop > start for each subject guarrantees that
    **    the 9 are a subset of the 27). 
    **  A basic rule is to remove subjects from the sums as soon as possible
    **    and add them as late as possible: (a-b) +c vs (a+c)-b when b might
    **    be large.
    **  Algorithm: 
    **    look ahead to find the next death time 'dtime'
    **    remove anyone no longer at risk from the sums, finish calculating resid
    **    add anyone newly at risk to the sums
    **    update the cumulative hazard
    **    everything resets at the end of a stratum.
    **  The sort order is from large time to small, so we encounter a subject's
    **    ending time first, then their start time.
    **  The martingale residual for a subject is 
    **     status - (cumhaz at tstop - cumhaz at tstart+0)*score
    **  If there is an observation that is never in a risk set, i.e., the
    **   (time1, time2) pair overlaps no events, the parent routine will have
    **   sorted them to the end, set nused < nr, and their initial residual
    **   of zero never is updated.
    */
    person1 =0;
    denom =0;
    cumhaz =0;
    istrat = strata[sort1[0]];
    for (person2=0; person2 <nused; ) {
	/* find the next event time */
	for (k=person2; k<nused; k++) {
	    p2 = sort2[k];
	    if (strata[p2] != istrat) {
		/* start of a new stratum */
		for (; person1< nused; person1++) {
		    p1 = sort1[person1];
		    if (strata[p1] != istrat) break;
		    resid[p1] -= cumhaz * score[p1];
		}
		cumhaz =0;
		denom = 0;
		istrat= strata[p2];
	    }
	    if (event[p2] >0) {
		dtime = tstop[p2];
		break;
	    }
	}
	if (k==nused) break;

	/* 
	** Remove those subjects whose start time is to the right
	**  from the risk set, and finish computation of their residual
	*/
	for (;  person1 <nused; person1++) {
	    p1 = sort1[person1];
	    if (tstart[p1] < dtime || strata[p1] != istrat) break;
	    if (atrisk[p1] == 1) {
		denom -= score[p1] * weight[p1];
		resid[p1] -= cumhaz * score[p1];
            }
	}
	
	/*
        **  Add in new subjects 
        */
	deaths =0;
	e_denom =0;
	wtsum =0;
	for (k=person2; k< nused; k++) {
	    p2 = sort2[k];
	    if (tstop[p2] < dtime || strata[p2] != istrat) break;

	    if (event[p2] ==1) {  /* stop[p2] will be = to dtime */
		atrisk[p2] =1;
		resid[p2] = 1.0 + cumhaz * score[p2];
		deaths ++;
		denom += score[p2] * weight[p2];
		e_denom += score[p2] * weight[p2];
		wtsum += weight[p2];
	    }
	    else if (tstart[p2] < dtime) { /* start > stop > dtime possible */
		denom += score[p2] * weight[p2];
		atrisk[p2] =1;
		resid[p2] = cumhaz * score[p2];
	    }
	}		    
	
	/* compute the increment in hazard 
	** hazard = usual increment
	** e_hazard = efron increment, for tied deaths only
	*/
	if (method==0 || deaths==1) { /* Breslow */
	    hazard = wtsum/denom;
	    person2 = k;
	}
	else { /* Efron */
	    hazard =0;
	    e_hazard =0;  /* hazard experienced by a tied death */
	    wtsum /= deaths;   
	    for (i=0; i <deaths; i++) {
		temp = i/deaths;
		hazard += wtsum/(denom - temp*e_denom);
		e_hazard += wtsum * (1-temp)/(denom - temp*e_denom);
	    }

	    /* deaths don't get the full hazard increment */
	    temp = hazard - e_hazard;  /* would be 0 for Breslow */
	    for (; person2 <k; person2++) {
		p2 = sort2[person2];
		if (event[p2] >0) resid[p2] += temp*score[p2];
	    }
	}
	cumhaz += hazard;
    }

    /* finish up the last few */
    for (; person1< nused; person1++) {
	p1 = sort1[person1];
	if (atrisk[p1]==1) resid[p1] -= cumhaz * score[p1];
    }

    UNPROTECT(1);
    return(resid2);
}
	
	    
	    
