/*
** Compute the martingale residual for a Cox model.
** This routine does the same work as coxmart, except
**    it expects data in inverse time order
**    only does the Breslow method
**    exists for the sake of coxexact.fit
**
** Input
**      n       number of subjects
**      time    vector of times
**      status  vector of status values
**      score   the vector of subject scores, i.e., exp(beta*z)
**      strata  is =1 for the first obs of a strata
**      wt      case weights
** Output
**      the residual for each subject
*/
#include "survS.h"
#include "survproto.h"
void coxmart2(Sint   *sn,     double *time, 
	     Sint   *status, Sint   * strata,   double *score, 
	     double *wt,     double *resid)
    {
    int i,j;
    int n;
    double deaths, denom;
    double expected, hazard;

    n = *sn;

    /*
    **	Accumulate the weighted score in reverse time order (data order)
    **  Temporarily save the resulting hazard in the residual vector
    */
    denom =0;
    for (i=0; i<n; ) {
	if (strata[i] ==1) {  /*first obs of a new strata */
	    denom=0;  
	    }

	deaths =0;  /* number of deaths at this time */
	/* 
	** for all obs at this time, add up the #deaths and weighted score 
	**  Watch out for the odd case that the last observation of one
	**  strata and the first obs of the next strata have the same value
	**  for time -- these don't count as tied times.
	*/
	denom += score[i]*wt[i];
	deaths += status[i]* wt[i];
	for (j=i+1; j<n && time[j]==time[i] && strata[j]==0; j++) {
	    denom += score[j]*wt[j];
	    deaths += status[j] * wt[j];
	    }
	hazard = deaths/denom; /* Breslow estimate of hazard */
	resid[j-1] = hazard;  

	i = j; /* forward to the next unique time */
	}

    /*
    ** pass 2: accumulate the hazard from smallest time to largest
    */
    expected =0;
    for (i= n-1; i>=0; i--) {
	expected += resid[i];
	resid[i] = status[i] - score[i]*expected;
	if (strata[i] ==1) expected=0;  /* last obs of a strata */
	}
    }
