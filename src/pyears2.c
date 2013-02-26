/*
**  Person-years calculations.
**     same as pyears1, but no expected rates
**
**  Input:
**      n       number of subjects
**      ny      number of columns of y.
**      doevent does y have an 'events' column?  1=yes, 0=no
**              if ny=2 and doevent=1, then "start" is missing.
**      y[3,n]  contains start, stop, and event for each subject
**      wt      contains the case weights
**
**  output table's description
**      odim        number of dimensions
**      ofac[odim]  1=is a factor, 0=continuous (time based)
**      odims[odim] the number of rows, columns, etc
**      ocut[]      for each non-factor dimension, the odim[i]+1 cutpoints
**                        that define the intervals; concatonated.
**      odata[odim, n]  the subject data-- where each indexes into the
**                        expected table, at time 0.
**
** Output:
**      pyears     output table of person years
**      pn         number of observations that contribute to each cell
**      pcount     number of events
**      offtable   total person years that did not fall into the output table
**
** Scratch     allocated on the fly
**      scratch[edim]
*/
#include "survS.h"
#include "survproto.h"

/* names that begin with "s" will be re-declared in the main body */
void pyears2(Sint   *sn,      Sint   *sny,   Sint   *sdoevent, 
	     double *sy,      double *wt,    Sint   *sodim,    Sint   *ofac, 
	     Sint   *odims,   double *socut, double *sodata,
	     double *pyears,  double *pn,    double *pcount, 
	     double *offtable)
    {

    int i,j;
    int     n,
	    ny,
	    doevent,
	    odim;
    double  *start,
	    *stop,
	    *event,
	    **ocut,
	    **odata;
    double  *data;
    double  timeleft,
	    thiscell;
    int     index;
    int     dostart;
    int     d1;    /* a dummy */
    double  d2;    /* a dummy for pystep */
    double  eps;   /* round off protection */

    n = *sn;
    ny= *sny;
    doevent = *sdoevent;
    odim = *sodim;
    start = sy;

    if (ny==3 || (ny==2 && doevent==0)) {
	/* each subject has a "start" time */
	stop = sy +n;
	dostart =1;
	}
    else   {
	/* followup starts at time 0 */
	stop  = sy;
	dostart =0;
	}
    event = stop +n;
    odata = dmatrix(sodata, n, odim);
    data  = (double *) ALLOC(odim, sizeof(double));
    /*
    ** will be a ragged array
    */
    ocut = (double **)ALLOC(odim, sizeof(double *));
    for (i=0; i<odim; i++) {
	ocut[i] = socut;
	if (ofac[i]==0) socut += odims[i] +1;
	}

    /*
    ** Set the round off error to min(time[time>0]) * 1e-8
    **   The events are counted in the last cell to which person years are
    **   added in the while() loop below.  We don't want to "spill over" into
    **   a next (incorrect) cell due to accumulated round off, in the case
    **   that a subjects fu time exactly matches one of the cell boundaries.
    */
    eps =0; /* guard against the rare case that all(time==0) */
    for (i=0; i<n; i++) {
	if (dostart==1) timeleft = stop[i] - start[i];
	else timeleft= stop[i];
	if (timeleft >0) {
	    eps = timeleft;  /* starting guess for min = first non-zero value*/
	    break;
	    }
	}
    for (; i<n; i++) {
	if (dostart==1) timeleft = stop[i] - start[i];
	else timeleft= stop[i];
	if ((timeleft >0) && (timeleft < eps)) eps = timeleft;
	}
    eps *= 1e-8;

    *offtable =0;
    for (i=0; i<n; i++) {
	/*
	** initialize
	** "data" will be the vector of starting values for each subject
	**    for a factor variable this is just the cell number (1,2,3,...)
	**    for a continuous one it is the value, which will be matched to
	**       the ocuts list, by pystep, to figure out a cell number
	*/
	for (j=0; j<odim; j++) {
	    if (ofac[j] ==1 || dostart==0) data[j] = odata[j][i];
	    else                           data[j] = odata[j][i] + start[i];
	    }
	if (dostart==1) timeleft = stop[i] - start[i];
	else            timeleft = stop[i];
	/*
	** add up p-yrs
	**   d1 and d2 are only relevant to rate tables, so ignored here
	** If there are only factor variables, the loop below will finish in
	**   one iteration.  Otherwise, the value "thiscell" is the amount of
	**   person-years in the current cell, up to the next cell boundary
	**   (or the amount of time left, whichever is lesser).
	*/
	if (timeleft <=eps && doevent) {
	    /* we have to call pystep at least once to set the index */
	    pystep(odim, &index, &d1, &d2, data, ofac, odims, ocut, 1, 0);
	    }

	while (timeleft > eps) {
	    thiscell = pystep(odim, &index, &d1, &d2, data, ofac, odims, ocut,
				    timeleft, 0);
	    if (index >=0) {
		pyears[index] += thiscell * wt[i];
		pn[index] += 1;
		}
	    else *offtable += thiscell * wt[i];

	    for (j=0; j<odim; j++)
		if (ofac[j] ==0) data[j] += thiscell;
	    timeleft -= thiscell;
	    }
	if (index >=0 && doevent) pcount[index] += event[i] * wt[i];
	}
    }
