/*
**  Person-years calculations, in its most general
**
**  Input:
**      n       number of subjects
**      ny      number of columns of y.
**      doevent does y have an 'events' column?  1=yes, 0=no
**              if ny=2 and doevent=1, then "start" is missing.
**      method  if =1 do expected number of events, else expected number
**                    of person years
**      y[3,n]  contains start, stop, and event for each subject
**      weight     case weights
**
**    expected table
**      edim        number of dimensions of the expected table
**      efac[edim]  1=is a factor, 0=continuous (time based)
**                    >=2  special handling for US "calendar year"
**      edims[edim] the number of rows, columns, etc
**      ecut[ ]     the starting points for each non-factor dimension,
**                          strung together.
**      expect      the actual table of expected rates
**      edata[edim, n]  the subject data-- where each indexes into the
**                        expected table, at time 0.
**
**   output table's description
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
**      pexpect    expected number of events
**      offtable   total person years that did not fall into the output table
**
** Scratch   -- allocated on the fly
**      scratch[edim + odim]
*/
#include<math.h>
#include "survS.h"
#include "survproto.h"

/* names that begin with "s" will be re-declared in the main body */
void pyears1(Sint   *sn,      Sint   *sny,      Sint   *sdoevent, 
	     double *sy,      double *weight,       
             Sint   *sedim,   Sint   *efac, 
	     Sint   *edims,   double *secut,    double *expect, 
	     double *sedata,  Sint   *sodim,    Sint   *ofac, 
	     Sint   *odims,   double *socut,    Sint   *smethod, 
	     double *sodata,  double *pyears,   double *pn, 
	     double *pcount,  double *pexpect,  double *offtable)
    {

    int i,j;
    int     n,
	    ny,
	    doevent,
            method,
	    edim,
	    odim;
    double  *start,
	    *stop,
	    *event,
	    **ecut,
	    **ocut,
	    **edata,
	    **odata;
    double  *data,
	    *data2;
    double  timeleft,
	    thiscell,
	    etime,
	    et2;
    int     index,
	    indx, indx2;
    double  lwt;   /*this variable is returned by pystep, and controls the
		   "on the fly" linear interpolation done for the calandar
		   year dimension of rate tables */
    int     dostart;
    double  hazard, cumhaz;
    double  temp, lambda;
    double  eps;   /* protection against accumulated round off */

    n = *sn;
    ny= *sny;
    doevent = *sdoevent;
    method  = *smethod;
    edim = *sedim;
    odim = *sodim;
    start = sy;
    if (ny==3 || (ny==2 && doevent==0)) {
	stop = sy +n;
	dostart =1;
	}
    else   {
	stop  = sy;
	dostart =0;
	}
    event = stop +n;
    edata = dmatrix(sedata, n, edim);
    odata = dmatrix(sodata, n, odim);
    i=edim + odim;
    data  = (double *) ALLOC(i, sizeof(double));
    data2 = data + odim;
    /*
    ** ecut and ocut will be ragged arrays
    */
    ecut = (double **)ALLOC(edim, sizeof(double *));
    for (i=0; i<edim; i++) {
	ecut[i] = secut;
	if (efac[i]==0)     secut += edims[i];
	else if(efac[i] >1) secut += 1 + (efac[i]-1)*edims[i];
	}

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
	    eps = timeleft;
	    break;
	    }
	}
    for (; i<n; i++) {
	if (dostart==1) timeleft = stop[i] - start[i];
	else timeleft= stop[i];
	if (timeleft >0 && timeleft < eps) eps = timeleft;
	}
    eps *= 1e-8;

    *offtable =0;
    for (i=0; i<n; i++) {
	/*
	** initialize
	*/
	for (j=0; j<odim; j++) {
	    if (ofac[j] ==1 || dostart==0) data[j] = odata[j][i];
	    else                           data[j] = odata[j][i] + start[i];
	    }
	for (j=0; j<edim; j++) {
	    if (efac[j] ==1 || dostart==0) data2[j] = edata[j][i];
	    else                           data2[j] = edata[j][i] + start[i];
	    }
	if (dostart==1) timeleft = stop[i] - start[i];
	else timeleft= stop[i];

	cumhaz=0;
	/*
	** add up p-yrs
	*/
	if (timeleft <=eps && doevent) {
	    /* we have to call pystep at least once to set the indices */
	    pystep(odim, &index, &indx2, &lwt, data, ofac, odims, ocut, 1, 0);
	    }

	while (timeleft > eps) {
	    thiscell = pystep(odim, &index, &indx2, &lwt, data, ofac, odims,
				     ocut, timeleft, 0);
	    if (index >=0) {
		pyears[index] += thiscell * weight[i];
		pn[index] += 1;

		/* expected calc */
		etime = thiscell;
		hazard=0;
		temp =0;
		while (etime >0) {
		    /*
		    ** The hazard or survival curve (temp) calculated within
		    ** this loop don't depend on the case weight --- the
		    ** whole loop is only for one person, and hazard is a
		    ** function of time alone.  Once computed, however, the
		    ** total hazard added into the expected table 
		    ** is weighted.
		    */
		    et2 = pystep(edim, &indx, &indx2, &lwt, data2, efac,
				 edims, ecut, etime, 1);
		    if (lwt <1) lambda = (lwt*expect[indx] +
						     (1-lwt)*expect[indx2]);
		    else       lambda =  expect[indx];
		    if (method==0)
			temp += exp(-hazard)*(1-exp(-lambda*et2))/ lambda;
		    hazard += lambda * et2;

		    for (j=0; j<edim; j++)
			if (efac[j] !=1) data2[j] += et2;
		    etime -= et2;
		    }
		if (method==1) pexpect[index] += hazard * weight[i];
		else           pexpect[index] += exp(-cumhaz)*temp * weight[i];
		cumhaz += hazard;
		}
	    else  {
		*offtable += thiscell * weight[i];
		for (j=0; j<edim; j++)
		    if (efac[j] !=1) data2[j] += thiscell;
		}

	    for (j=0; j<odim; j++)
		if (ofac[j] ==0) data[j] += thiscell;
	    timeleft -=thiscell;
	    }
	if (index >=0 && doevent) pcount[index] += event[i] * weight[i];
	}
    }
