/*
** The fastest possible calculation (at least the fastest I've come up with)
** for martingale residuals from a (start, stop] model.  
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

SEXP agmart3(SEXP nused2,  SEXP surv2,  SEXP score2, SEXP weight2, 
	     SEXP strata2, SEXP sort12, SEXP sort22, SEXP method2) {

    int i, k;
    int p, p1, istrat, indx2;
    double deaths, denom, e_denom;
    double hazard, e_hazard, cumhaz;
    double temp, dtime;
    double wtsum;
    int nr, person, nused;
 
    /* pointers to the input data */
    double *start, *stop, *event;
    double *weight, *score;
    int    *sort1, *sort2, *strata;

    int method;  /* integer version of input */

    /* output */
    SEXP resid2;
    double *resid;

    nused = asInteger(nused2);
    nr = nrows(surv2);
    method = asInteger(method2);
    start = REAL(surv2);
    stop  = start +nr;
    event = stop +nr;
    weight= REAL(weight2);
    score = REAL(score2);
    sort1 = INTEGER(sort12);
    sort2 = INTEGER(sort22);
    strata= INTEGER(strata2);

    PROTECT(resid2 = allocVector(REALSXP, nr));
    resid = REAL(resid2);
    for (k=0; k<nr; k++) resid[k] =0;

    /*
    **  'person' walks through the the data from 1 to nused,
    **     sort1[0] points to the largest stop time, sort1[1] the next, ...
    **  'time' is a scratch variable holding the time of current interest
    **  'indx2' walks through the start times.  It will be smaller than 
    **    'person': if person=27 that means that 27 subjects have stop >=time,
    **    and are thus potential members of the risk set.  If 'indx2' =9,
    **    that means that 9 subjects have start >=time and thus are NOT part
    **    of the risk set.  (stop > start for each subject guarrantees that
    **    the 9 are a subset of the 27). 
    **  A basic rule is to remove subjects from the sums as soon as possible
    **    and add them as late as possible: (a-b) +c vs (a+c)-b when b might
    **    be large.
    **  Algorithm: 
    **    look ahead to find the next death time 'dtime'
    **    finish the calc for anyone no longer at risk at that time, then remove
    **    add anyone newly at risk to the sums
    **    update the cumulative hazard
    **    everything resets at the end of a stratum.
    **  The sort order is from large time to small, so we encounter a subject's
    **    ending time first, then their start time.
    **  The martingale residual for a subject is 
    **     status - (cumhaz at end of their interval - cumhaz at start)*score
    **  If there is an observation that is never in a risk set, i.e., the
    **   (time1, time2) pair overlaps no events, the parent routine will have
    **   sorted them to the end, set nused < nr, and their initial residual
    **   of zero never is updated.
    */
    indx2 =0;
    denom =0;
    cumhaz =0;
    istrat = strata[sort1[person]];
    for (person=0; person <nused; ) {
	/* find the next event time */
	for (k=person; k<nused; k++) {
	    p = sort2[k];
	    if (strata[p] != istrat) {
		/* start of a new stratum */
		for (; indx2< nused; indx2++) {
		    p1 = sort1[indx2];
		    if (strata[p1] != istrat) break;
		    resid[p1] -= cumhaz * score[p1];
		}
		cumhaz =0;
		denom = 0;
		istrat= strata[p];
	    }
	    if (event[p] >0) break;
	}

	if (k==nused) {
	    /* no event found, all done */
	    for (; indx2< nused; indx2++) {
		p1 = sort1[person];
		resid[p1] -= cumhaz * score[p1];
	    }
	    person = nused;
	    break;
	}	
	else dtime = stop[p];

	/* 
	** Remove those subjects whose start time is to the right
	**  from the risk set, and finish computation of their residual
	*/
	for (;  indx2 <nused; indx2++) {
	    p1 = sort1[indx2];
	    if (start[p1] <= dtime || strata[p1] != istrat) break;
	    denom -= score[p] * weight[p1];
	    resid[p] -= cumhaz * score[p1];
	}
	
	/*
        **  Add in new subjects 
        */
	deaths =0;
	e_denom =0;
	for (k=person; k< nused; k++) {
	    p = sort2[k];
	    if (stop[p] < dtime || strata[p] != istrat) break;

	    denom += score[p] * weight[p];
	    resid[p] = cumhaz * score[p];
	    if (event[p] ==1) {
		deaths ++;
		e_denom += score[p] * weight[p];
		wtsum += weight[p];
	    }
	}		    
	
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
	    for (i=0; i <deaths; i++) {
		temp = i/deaths;
		hazard += wtsum/(denom - temp*e_denom);
		e_hazard += wtsum * (1-temp)/(denom - temp*e_denom);
	    }
	}
	cumhaz += hazard;

	/*
	** Update the hazard for deaths
        */
	temp = hazard - e_hazard;  /* will be 0 for Breslow */
	for (; person <k; person++) {
	    p = sort2[person];
	    if (event[p] >0) resid[p] += 1 + temp*score[p];
	}
    }

    UNPROTECT(1);
    return(resid2);
}
	
	    
	    
