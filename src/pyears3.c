/*
**  Person-years calculations, leading to expected survival for a cohort.
**    The output table depends only on factors, not on continuous.
**
**  Input:
**      death        1=conditional surv, 0=cohort
**      n            number of subjects
**
**    expected table
**      edim        number of dimensions of the expected table
**      efac[edim]  1=is a factor, 0=continuous (time based)
**      edims[edim] the number of rows, columns, etc
**      ecut[sum(edims)]  the starting point (label) for each dimension.
**                          if it is a factor dim, will be 1:edims[i]
**      expect      the actual table of expected rates
**
**    subject data
**      group
**      x[edim, n]  the first column is the group, the rest are where each
**                      subject indexes into the expected table, at time 0
**      y[n]         the time at risk for each subject
**
**    control over output
**      ntime           the number of time points desired
**      ngrp            the number of patient groups
**      times           the list of output times
**
**    Output
**      esurv[ntime,ngrp]   conditional survival
**      nsurv[ntime,ngrp]   number of subjects per cell of "esurv"
*/
#include <math.h>
#include "survS.h"
#include "survproto.h"

/* names that begin with "s" will be re-declared in the main body */
void pyears3(Sint   *sdeath,    Sint   *sn,    Sint   *sedim, 
	     Sint   *efac,      Sint   *edims, double *secut, 
	     double *expect,    Sint   *grpx, double *sx,    double *y, 
	     Sint   *sntime,    Sint   *sngrp, double *times,
	     double *esurv,     Sint   *nsurv)
    {

    int i,j,k;
    int     n,
	    death,
	    edim,
	    ngrp,
	    ntime;
    double  **x;
    double  *data2;
    double  **ecut;
    double  hazard,   /*cum hazard over an interval */
	    cumhaz;   /*total hazard to date for the subject */
    double  timeleft,
	    thiscell,
	    etime,
	    et2;
    int     index,
	    indx,
	    indx2;
    double  wt;
    double  *wvec;    /* vector of weights needed for unconditional surv */
    int     group;
    double  time;

    death = *sdeath;
    n = *sn;
    edim = *sedim;
    ntime = *sntime;
    ngrp  = *sngrp;
    x     = dmatrix(sx, n, edim);
    data2 = (double *)ALLOC(edim+1, sizeof(double));
    wvec  = (double *)ALLOC(ntime*ngrp, sizeof(double));
    for (j=0; j<ntime*ngrp; j++) wvec[j] =0;

    /*
    ** ecut will be a ragged array
    */
    ecut = (double **)ALLOC(edim, sizeof(double *));
    for (i=0; i<edim; i++) {
	ecut[i] = secut;
	if (efac[i]==0)     secut += edims[i];
	else if(efac[i] >1) secut += 1 + (efac[i]-1)*edims[i];
	}

    for (i=0; i<n; i++) {
	R_CheckUserInterrupt();  /*check for control-C */
	/*
	** initialize
	*/
	cumhaz =0;
	for (j=0; j<edim; j++) data2[j] = x[j][i];
	timeleft = y[i];
	group = grpx[i] -1;
	time =0;      /*change this later to an input paramter, i.e., start */

	/*
	** add up hazard
	*/
	for (j=0; j<ntime && timeleft >0; j++) {
	    thiscell = times[j] - time;
	    if (thiscell > timeleft) thiscell = timeleft;
	    index =j + ntime*group;

	    /* expected calc 
	    **  The wt parameter only comes into play for older style US rate
	    **   tables, where pystep does interpolation.
	    ** Each call to pystep moves up to the next 'boundary' in the
	    **  expected table, data2 contains our current position therein
	    */
	    etime = thiscell;
	    hazard =0;
	    while (etime >0) {
		et2 = pystep(edim, &indx, &indx2, &wt, data2, efac,
			     edims, ecut, etime, 1);
		if (wt <1) hazard+= et2*(wt*expect[indx] +(1-wt)*expect[indx2]);
		else       hazard+= et2* expect[indx];
		for (k=0; k<edim; k++)
		    if (efac[k] !=1) data2[k] += et2;
		etime -= et2;
/*
printf("indx=%d, time=%5.1f, rate1=%6e, rate2=%6e, wt=%3.1f\n", 
       indx, et2, expect[indx], expect[indx2], wt);
*/
		}
/*printf("index=%d, hazard=%6e, cumhaz=%6e\n", index, hazard, cumhaz); */
	    if (times[j]==0) {
		wvec[index]=1;
		if (death==0) esurv[index]=1;
		else          esurv[index]=0;
		}
	    else if (death==0) {
		esurv[index] += exp(-(cumhaz+hazard)) * thiscell;
		wvec[index]  += exp(-cumhaz) * thiscell;
		}
	    else {
		esurv[index] += hazard * thiscell;
		wvec[index] +=  thiscell;
		}
	    nsurv[index] ++;
	    cumhaz += hazard;

	    time  += thiscell;
	    timeleft -= thiscell;
	    }
	}

    for (i=0; i<ntime*ngrp; i++) {
/*
printf("i=%3d, esurv=%6e, wvec=%6e, death=%d\n", i, esurv[i], wvec[i], death);
*/
	if (wvec[i]>0) {
	    if (death==0) esurv[i] /= wvec[i];
	    else          esurv[i] = exp(-esurv[i]/wvec[i]);
	    }
	else if (death!=0) esurv[i] = exp(-esurv[i]);
	}
    }
