/*
** The fastest possible calculation (at least the fastest I've come up with)
** for martingale residuals from a (start, stop] model.  Motivated by the need
** of coxme to do this multiple times.
**
** Input
**      surv    start, stop, event matrix of survival data
**      score   the vector of subject scores, i.e., exp(beta*z)
**      weight  case weights
**      strata   sizes of the strata
**      sort     sort order for the obs: col1 = stop times, 2 = start times
**         time within strata for both
**      method  will be ==1 for the Efron method
** Output
**      martingale residual
*/
#include <math.h>
#include "survS.h"
#include "survproto.h"

SEXP agmart3(SEXP surv2, SEXP score2, SEXP weight2, SEXP strata2,
	     SEXP sortx, SEXP method2) {

    int k, ksave;
    int p, istrat, indx2;
    double deaths, denom, e_denom;
    double hazard, e_hazard, cumhaz;
    double temp, time;
    double wtsum;
    int n, person;
 
    /* pointers to the input data */
    double *start, *stop, *event;
    double *weight, *score;
    int    *sort1, *sort2, *strata;

    int method;  /* integer version of input */

    /* output */
    SEXP resid2;
    double *resid;

    n = nrows(surv2);
    method = asInteger(method2);
    start = REAL(surv2);
    stop  = start +n;
    event = stop +n;
    weight= REAL(weight2);
    score = REAL(score2);
    sort1 = INTEGER(sortx);
    sort2 = sort1 + n;
    strata= INTEGER(strata2);

    PROTECT(resid2 = allocVector(REALSXP, n));
    resid = REAL(resid2);

    /*
    **  'person' walks through the the data from 1 to n,
    **     sort1[0] points to the largest stop time, sort1[1] the next, ...
    **  'time' is a scratch variable holding the time of current interest
    **  'indx2' walks through the start times.  It will be smaller than 
    **    'person': if person=27 that means that 27 subjects have stop >=time,
    **    and are thus potential members of the risk set.  If 'indx2' =9,
    **    that means that 9 subjects have start >=time and thus are NOT part
    **    of the risk set.  (stop > start for each subject guarrantees that
    **    the 9 are a subset of the 27). 
    **  Basic algorithm: move 'person' forward, adding the new subject into
    **    the risk set.  If this is a new, unique death time, take selected
    **    old obs out of the sums, add in obs tied at this time, then update
    **    the cumulative hazard. Everything resets at the end of a stratum.
    **  The sort order is from large time to small, so we encounter a subject's
    **    ending time first, then their start time.
    **  The martingale residual for a subject is 
    **     status - (cumhaz at end of their interval - cumhaz at start)*score
    */
    istrat=0;
    indx2 =0;
    denom =0;
    cumhaz =0;
    for (person=0; person <n; ) {
	p = sort1[person];
	if (event[p] ==0) { /* censored */
	    denom += score[p] * weight[p];
	    resid[p] = cumhaz * score[p];
	    person++;
	} else {
	    time = stop[p]; /* found a new, unique death time */
	    /* 
	    ** Remove those subjects whose start time is to the right
	    **  from the risk set, and finish computation of their residual
	    */
	    for (;  indx2 <strata[istrat]; indx2++) {
		p = sort2[indx2];
		if (start[p] < time) break;
		denom -= score[p] * weight[p];
		resid[p] -= cumhaz * score[p];
	    }

	    /*
	    **	Add up over this death time, for all subjects
	    */
	    deaths =0;
	    e_denom =0;
	    wtsum =0;
	    for (k=person; k<strata[istrat]; k++) {
		p = sort1[k];
		if (stop[p]  < time) break;  /* only tied times */ 
		denom += score[p] * weight[p];
		if (event[p] ==1) {
		    deaths ++;
		    e_denom += score[p] * weight[p];
		    wtsum += weight[p];
		}
	    }
	    ksave = k;
	    
	    /* compute the increment in hazard 
	    ** hazard = usual increment
	    ** e_hazard = efron increment, for tied deaths only
	    */
	    if (method==0 || deaths==1) { /* Breslow */
		hazard = wtsum/denom;
		e_hazard = hazard;
	    }
	    else { /* Efron */
		hazard =0;
		e_hazard =0;  /* hazard experienced by a tied death */
		wtsum /= deaths;   
		for (k=0; k <deaths; k++) {
		    temp = k/deaths;
		    hazard += wtsum/(denom - temp*e_denom);
		    e_hazard += wtsum * (1-temp)/(denom - temp*e_denom);
		}
	    }

	    /* Give initial value to all intervals ending at this time
            ** If tied censors are sorted before deaths (which at least some
	    **  callers of this routine do), then the else below will never
	    **  occur.
            */
	    temp = cumhaz + (hazard -e_hazard);
	    for (; person < ksave; person++) {
		p = sort1[person];
		if (event[p] ==1) resid[p] = 1 + temp*score[p];
		else resid[p] = cumhaz * score[p];
		}
	    cumhaz += hazard;
	}

	/* clean up at the end of a strata */
	if (person == strata[istrat]) {
	    for (; indx2<strata[istrat]; indx2++) {
		p = sort2[indx2];
		resid[p] -= cumhaz * score[p];
	    }
	    cumhaz =0;
	    denom = 0;
	    istrat++;
	}
    }
    UNPROTECT(1);
    return(resid2);
}
	
	    
	    
