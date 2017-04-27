/* $Id: survdiff2.c 11166 2008-11-24 22:10:34Z therneau $ */
#include <math.h>
#include "survS.h"
#include "survproto.h"

void survdiff2(Sint   *nn,     Sint   *nngroup,    Sint   *nstrat, 
	       double *rho,    double *time,       Sint   *status, 
	       Sint   *group,  Sint   *strata,	   double *obs, 
	       double *exp,    double *var,        double *risk, 
	       double *kaplan)
    {
    register int i,j,k;
    int kk;
    int n, ngroup, ntot;
    int istart, koff;
    double km, nrisk, wt, tmp;
    double deaths;

    ntot = *nn;
    ngroup = *nngroup;
    istart=0; koff=0;
    for (i=0; i< ngroup*ngroup; i++)  var[i]=0;
    for (i=0; i< *nstrat*ngroup; i++) {
	obs[i]=0;
	exp[i]=0;
	}

    while (istart < ntot) {  /* loop over the strata */
	for (i=0; i<ngroup; i++) risk[i]=0;

	/* last obs of this strata */
	for (i=istart; i<ntot; i++)
	    if (strata[i]==1) break;
	n = i+1;


	/*
	** Compute the k-m, which is only needed if rho!=0
	**   We want it set up as a left-continuous function (unusual)
	*/
	if (*rho !=0){
	    km =1;
	    for (i=istart; i<n; ) {
		kaplan[i] = km;
		nrisk = n-i;
		deaths =status[i];
		for (j=i+1; j<n && time[j]==time[i]; j++) {
		    kaplan[j] = km;
		    deaths += status[j];
		    }
		km = km * (nrisk-deaths)/nrisk;
		i=j;
		}
	    }

	/*
	** Now for the actual test
	*/
	for (i=n-1; i>=istart; i--) {
	    if (*rho ==0) wt=1;
	    else          wt= pow(kaplan[i], *rho);

	    deaths = 0;
	    for (j=i; j>=istart && time[j]==time[i]; j--) {
		k = group[j]-1;
		deaths += status[j];
		risk[k] += 1;
		obs[k + koff] += status[j] *wt;
		}
	    i= j +1;
	    nrisk = n-i;

	    if (deaths>0) {  /* a death time */
		for (k=0; k<ngroup; k++)
		    exp[k+koff] += wt* deaths * risk[k] / nrisk;

		if (nrisk==1) continue;  /*only 1 subject, so no variance */
		kk =0;
		wt = wt*wt;
		for (j=0; j<ngroup; j++) {
		    tmp = wt* deaths* risk[j]* (nrisk-deaths)/(nrisk *(nrisk-1));
		    var[kk+j] += tmp;
		    for (k=0; k<ngroup; k++) {
			var[kk] -= tmp * risk[k] / nrisk;
			kk++ ;
			}
		    }
		}
	    }
	istart = n;
	koff += ngroup;
	}
    }
