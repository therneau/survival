/* $Id: coxscho.c 11376 2009-12-14 22:53:57Z therneau $
** Return the Schoenfeld residuals.
**
**  the input parameters are
**
**       nused        :number of people
**       nvar         :number of covariates
**       y(3,n)       :start, stop, and event for each subject
**       covar(nv,n)  :covariates for person i.
**                        Note that S sends this in column major order.
**       strata(n)    :marks the strata.  Will be 1 if this person is the
**                       last one in a strata.  If there are no strata, the
**                       vector can be identically zero, since the nth person's
**                       value is always assumed to be = to 1.
**       score(n)     :the risk score for the subject
**       method       : =1 if the Efron method was used
**
**  returned parameters
**       covar        :for each death, the row is replaced with the residuals
**
**  work arrays
**       a(nvar)
**       a2(nvar)
**       mean(nvar)
**
**  the 3 vectors a, a2 and mean are passed as a single
**    vector of storage, and then broken out.
**
**  the data must be sorted by ascending time within strata, deaths before
**          living within tied times.
*/
#include <math.h>
#include "survS.h"
#include "survproto.h"

void coxscho(Sint   *nusedx,    Sint   *nvarx,    double *y, 
	     double *covar2,    double *score,    Sint   *strata,  
	     Sint   *method2,   double *work)
{
    int i,k,person;
    int     nused, nvar;
    double **covar;
    double *a;
    double *a2;
    double *mean;
    double  denom, weight;
    double  time;
    double  temp;
    double     method;
    double   deaths;
    double efron_wt;
    double  *start,
	    *stop,
	    *event;

    nused = *nusedx;
    nvar  = *nvarx;
    method= *method2;
    /*
    **  Set up the ragged arrays
    */
    covar= dmatrix(covar2, nused, nvar);
    a = work;
    a2= a+nvar;
    mean = a2+nvar;
    start =y;
    stop  =y + nused;
    event =y + nused +nused;

    /*
    ** Now walk through the data
    */
    for (person=0; person<nused;) {
	if (event[person]==0) person++;
	else {
	    /*
	    ** compute the mean over the risk set and over the deaths (a & a2)
	    */
	    denom =0;
	    efron_wt =0;
	    for (i=0; i<nvar; i++) {
		a[i] =0;
		a2[i]=0;
		}
	    time = stop[person];
	    deaths=0;
	    for (k=person; k<nused; k++) {
		if (start[k] < time) {
		    weight = score[k];
		    denom += weight;
		    for (i=0; i<nvar; i++) {
			a[i] += weight*covar[i][k];
			}
		    if (stop[k]==time && event[k]==1) {
			deaths += 1;
			efron_wt += weight*event[k];
			for (i=0; i<nvar; i++) a2[i]+= weight*covar[i][k];
			}
		     }
		if (strata[k]==1) break;
		}

	    /*
	    ** Compute the mean at this time point
	    */
	    for (i=0; i<nvar; i++) mean[i] =0;
	    for (k=0; k<deaths; k++) {
		temp = method *k/deaths;
		for (i=0; i<nvar; i++)
		    mean[i] += (a[i] - temp*a2[i])/(deaths*(denom -temp*efron_wt));
		}
	    /*
	    ** Compute the residual(s) for this time point
	    */
	    for (k=person; k<nused && stop[k]==time; k++) {
		if (event[k]==1) {
		    for (i=0; i<nvar; i++) {
			covar[i][k] -= mean[i];
			}
		    }
		person++;
		if (strata[k]==1) break;
		}
	    }
	}
    }
