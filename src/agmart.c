/*
** Compute the martingale residual for a counting process model
**  This routine is only used by the exact partial likelihood routine
**  agexact.fit, which has such low usage that it and components have
**  not been coverted to faster algorithms as other routines have.
**
** Input
**      n       number of subjects
**      method  will be ==1 for the Efron method
**      start
**      stop    vector of (start, stop] times for the subjects
**      event   vector of status values
**      score   the vector of subject scores, i.e., exp(beta*z)
**      weight  case weights
**      strata  is =1 for the last obs of a strata
**
** Output
**      resid   martingale residual
**
** The martingale residual is more of a nuisance for the Efron method
**
**  Data must be sorted by time within strata, with events first
*/
#include <stdio.h>
#include "survS.h"
#include "survproto.h"

void agmart(Sint   *n,     Sint   *method,  double *start,   double *stop, 
	    Sint   *event, double *score,   double *wt,      Sint   *strata, 
	    double *resid)
    {
    int i,k;
    double deaths, denom, e_denom;
    double hazard, e_hazard;
    double temp, time;
    double wtsum;
    int nused;
    int person;

    nused = *n;
    strata[nused-1] =1;  /* Failsafe */

    for (i=0; i<nused; i++)  resid[i]=event[i];
    for (person=0; person<nused;) {
	if (event[person]==0) person++;
	else {
	    denom =0;
	    e_denom =0;
	    wtsum =0;
	    time = stop[person];
	    deaths=0;
	    for (k=person; k<nused; k++) {
		if (start[k] < time) {
		    denom += score[k]*wt[k];
		    if (stop[k]==time && event[k]==1) {
			deaths++;
			wtsum += wt[k];
			e_denom += score[k]*wt[k];
			}
		     }
		if (strata[k]==1) break;
		}

	    /*
	    ** Do "expected" for the risk set
	    */
	    hazard =0;
	    e_hazard=0;
	    wtsum /=deaths;
	    for (k=0; k<deaths; k++) {
		temp = *method *(k/deaths);
		hazard += wtsum/(denom - temp*e_denom);
		e_hazard += wtsum*(1-temp)/(denom - temp*e_denom);
		}
	    for (k=person; k<nused; k++) {
		if (start[k] < time) {
		    if (stop[k]==time && event[k]==1)
			    resid[k] -= score[k]*e_hazard;
		    else    resid[k] -= score[k]*hazard;
		    }
		if (stop[k]==time) person++;
		if (strata[k]==1) break;
		}
	    }
	}
    }
