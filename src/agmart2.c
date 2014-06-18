/*
** Compute the martingale residual for a counting process model
**   This routine has been superseded by agmart3
** Input
**      n       number of subjects
**      method  will be ==1 for the Efron method
**      start
**      stop    vector of (start, stop] times for the subjects
**      event   vector of status values
**      score   the vector of subject scores, i.e., exp(beta*z)
**      weight  case weights
**      nstrat       :number of strata
**      strata(nstrat):sizes of the strata
**      sort1        :sort order for the obs, using stop time, last to first
**      sort2        :sort order for the obs, using start time, last to first
**	haz	     :scratch vector of length 2 * #events
**		      actually, 2* max_unique_event_times_in_any_1_strata would
**		      suffice, but computing that bound isn't worth the bother
** Output
**      resid   martingale residual
**
** The martingale residual is more of a nuisance for the Efron method
*/
#include <stdio.h>
#include "survS.h"
#include "survproto.h"

void agmart2(Sint   *n,     Sint   *method,  double *start,   double *stop, 
	    Sint   *event,  Sint   *nstrat,  Sint *strata,    Sint *sort1,
	    Sint   *sort2,  double *score,   double *wt,      
	    double *resid,  double *haz)
    {
    int i, j, k, ksave;
    int p, istrat, indx2;
    double deaths, denom, e_denom;
    double hazard, e_hazard;
    double temp, time;
    double wtsum, *dtimes;
    int nused, ndeath;
    int person;
    int strata_start;

    nused = *n;
    j=0;
    for (i=0; i<nused; i++) {
	resid[i]=event[i];
	j += event[i];
	}
    dtimes = haz +j;    /*haz & dtimes were passed as one scratch vector */
    indx2 =0;
    denom =0;
    istrat=0;
    ndeath = 0;
    strata_start =0;
    for (person=0; person<nused;) {
	p = sort1[person];
	if (event[p]==0) {
	    denom += score[p]*wt[p];
	    person++;
	    }
	else {
	    e_denom =0;
	    wtsum =0;
	    time = stop[p];
	    deaths=0;			/* # tied deaths at this time */
	    for (k=person; k<strata[istrat]; k++) {
		p = sort1[k];
		if (stop[p] < time) break;
		denom += score[p] * wt[p];
		if (event[p]==1) {
		    deaths += 1;
		    e_denom += score[p]*wt[p];
		    wtsum += wt[p];
		    }
		}
	    ksave = k;

	    /*
	    ** subtract out the subjects whose start time is to the right
	    */
	    for (; indx2<strata[istrat]; indx2++) {
		p = sort2[indx2];
		if (start[p] < time) break;
		denom -= score[p] * wt[p];
		}

	    /*
	    ** At this point denom, e_denom, etc are updated.  They are
	    **   used to create the increment to the hazard:
	    **   hazard = usual increment
	    **   e_hazard = efron increment, for tied deaths only
	    */
	    hazard =0;
	    e_hazard=0;
	    wtsum /=deaths;
	    for (k=0; k<deaths; k++) {
		temp = *method *(k/deaths);
		hazard += wtsum/(denom - temp*e_denom);
		e_hazard += wtsum*(1-temp)/(denom - temp*e_denom);
		}
	    dtimes[ndeath] = time;	/* remember the death times */
	    haz[ndeath++] = hazard;

	    /*
	    ** Add this hazard increment to all whose intervals end just now
	    */
	    for (k=person-1; k>=strata_start; k--) { /*non-deaths */
		p = sort1[k];
		if (stop[p] > time) break;
		resid[p] -= score[p]*hazard;
		}
	    for (; person<ksave; person++) { /* deaths */
		p = sort1[person];
		resid[p] -= score[p]*e_hazard;
		}
	    }
	    
	if (person==strata[istrat]) {
	    /*
	    ** Now, walk through the list of subjects, adding in the
	    **  accrued hazard at non-death times (times strictly inside
	    **  each subject's interval).
	    */
	    k=0;
	    for (i=strata_start; i<person; i++) {
		p = sort1[i];
		for (; k< ndeath && dtimes[k] >= stop[p]; k++);  
		for (j=k; j<ndeath; j++)
		    if (start[p] < dtimes[j]) resid[p] -= score[p]*haz[j];
		}
	    istrat++;
	    denom=0;
	    ndeath=0;
	    indx2 =  person;
	    strata_start=person;
	    }
	}
    }
