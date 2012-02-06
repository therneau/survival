/*
**  Person-years calculations, leading to expected survival for a cohort.
**    The output table depends only on factors, not on continuous.
**    This version converted to .Call syntax for memory savings
**
**  Input:
**      death        1=conditional surv, 0=cohort
**
**    expected table, a multi-way array
**      efac[edim]  1=is a factor, 0=continuous (time based)
**      edims[edim] the dimension vector of the table; edim is its length
**      ecut[sum(edims)]  the starting point (label) for each dimension.
**                          if it is a factor dim, will be 1:edims[i]
**      expect      the actual table of expected rates
**
**    subject data
**      grpx[n]     how patients are grouped into curves
**      x[edim, n]  where each subject indexes into the expected table 
**                       at time 0, n= number of subjects
**      y[n]         the time at risk for each subject
**
**    control over output
**      ngrp  number of output groups (the max value of grpx)
**      times[ntime]    the list of output times
**
**    Output
**      esurv[ntime,ngrp]   conditional survival
**      nsurv[ntime,ngrp]   number of subjects per cell of "esurv"
*/
#include <math.h>
#include "survS.h"
#include "survproto.h"

/* my habit is to name a S object "charlie2" and the pointer
**  to the contents of the object "charlie"; the latter is
**  used in the computations
*/
SEXP pyears3b(SEXP   death2,    SEXP   efac2,   SEXP edims2,
	      SEXP   ecut2,     SEXP   expect2, SEXP grpx2,
	      SEXP   x2, 	SEXP   y2,      SEXP times2,
	      SEXP   ngrp2) {
    int i,j,k;
    int     n,
	    death,
	    edim,
	    ngrp,
	    ntime;
    double  **x;
    double  *data2;
    double  **ecut, *etemp;
    double  hazard,   /*cum hazard over an interval */
	    cumhaz;   /*total hazard to date for the subject */
    double  timeleft,
	    thiscell,
	    etime,
	    time,
	    et2;
    int     index,
	    indx,
	    indx2;
    double  wt;
    double  *wvec;    /* vector of weights needed for unconditional surv */
    int     group;

    int	    *efac, *edims, *grpx;
    double  *expect, *y, *times;
    SEXP    esurv2, nsurv2, rlist, rlistnames;
    double  *esurv;
    int     *nsurv;
    

    /* 
    ** copies of input arguments
    */
    death = asInteger(death2);
    ngrp  = asInteger(ngrp2);
    efac  = INTEGER(efac2);
    edims = INTEGER(edims2);
    edim  = LENGTH(edims2);
    expect= REAL(expect2);
    grpx  = INTEGER(grpx2);
    
    n     = LENGTH(y2);
    x     = dmatrix(REAL(x2), n, edim);
    y     = REAL(y2);
    times = REAL(times2);
    ntime = LENGTH(times2);

    /* scratch space */
    data2 = (double *)ALLOC(edim+1, sizeof(double));
    wvec  = (double *)ALLOC(ntime*ngrp, sizeof(double));
    for (j=0; j<ntime*ngrp; j++) wvec[j] =0;

    /*
    ** Set up ecut index as a ragged array
    */
    ecut = (double **)ALLOC(edim, sizeof(double *));
    etemp = REAL(ecut2);
    for (i=0; i<edim; i++) {
	ecut[i] = etemp;
	if (efac[i]==0)     etemp += edims[i];
	else if(efac[i] >1) etemp += 1 + (efac[i]-1)*edims[i];
	}

    /* 
    ** Create output arrays
    */
    PROTECT(esurv2 = allocVector(REALSXP, ntime*ngrp));
    esurv = REAL(esurv2);
    PROTECT(nsurv2 = allocVector(INTSXP, ntime*ngrp));
    nsurv = INTEGER(nsurv2);
    for (i=0; i<(ntime*ngrp); i++) {
	esurv[i] =0.;
	nsurv[i] =0;
	}

    /* compute */
    for (i=0; i<n; i++) {
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
    
    /*
    ** package the output
    */
    PROTECT(rlist = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(rlist,0, esurv2);
    SET_VECTOR_ELT(rlist,1, nsurv2);
    
    PROTECT(rlistnames= allocVector(STRSXP, 2));
    SET_STRING_ELT(rlistnames, 0, mkChar("surv"));
    SET_STRING_ELT(rlistnames, 1, mkChar("n"));
    setAttrib(rlist, R_NamesSymbol, rlistnames);

    unprotect(4);
    return(rlist);
    }
