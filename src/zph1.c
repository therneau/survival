/* 
** a 1 step variant of coxph5 for the zph test.
** gt = a function of time, with one element per unique death time
**   (within strata)
**  only compute the information needed for the score test,
**  which is the vector \sum g(t)* (x_i - xbar) of p elements and the 
**  full 2p by 2p information matrix
**
** gt = g(time) vector, centered
** y  = survival
** covar = X matrix of covariates
** eta   = risk score X beta
** weights = case weights
** strata  = stratum, integers of 1, 2, ...
** method = 2 for Efron, 1 for Breslow
** sort   = ordering vector
*/

#include <math.h>
#include "survS.h"
#include "survproto.h"

SEXP zph1(SEXP gt2,    SEXP y2, 
	  SEXP covar2, SEXP eta2,  SEXP weights2,
	  SEXP strata2,SEXP method2, SEXP sort2) {

    int i,j, k, kk, person, ip;
    int cstrat;   /* current stratum*/
    double temp, temp2, tmean;
    double *dtemp, timewt=1;
    double *a, *a2, **cmat, **cmat2;
    double denom=0, dtime, ndead, denom2;
    double risk, deadwt, wtave;
    int nprotect;

    /* scalar input arguments and counts*/
    int  nused, nvar, nevent, nstrat;
    int  method;

    /* input vectors */
    double *gt, *eta, *weights, *status;
    int    *strata, *sort;
    double **covar, *xtime;

    /* returned objects */
    SEXP rlist;
    double *u, **imat, **schoen;
    int **used;
    static const char *rnames[]={"u", "imat", "schoen", "used", ""};

    /* get local copies of input args */
    nused = nrows(y2);
    nvar  = ncols(covar2);
    method = asInteger(method2);

    xtime =   REAL(y2);
    status =  xtime + nused;
    weights = REAL(weights2);
    strata = INTEGER(strata2);
    gt =     REAL(gt2);
    eta =    REAL(eta2);
    sort=    INTEGER(sort2);

    /* 
    ** count up the number of events and the number of strata
    **  strata are numbered from 0
    */
    nevent =0;
    nstrat =1;
    for (i=0; i< nused; i++) {
	nevent += status[i];
	if (strata[i] >= nstrat) nstrat= strata[i]+1;
    }

    /*
    **  Set up the ragged arrays and scratch space
    */
    nprotect =0;
    if (MAYBE_REFERENCED(covar2)){
	nprotect++;
	covar = dmatrix(REAL(PROTECT(duplicate(covar2))), nused, nvar);
	}
    else covar = dmatrix(REAL(covar2), nused, nvar);

    nprotect++;
    PROTECT(rlist = mkNamed(VECSXP, rnames));
    u = REAL(SET_VECTOR_ELT(rlist, 0, allocVector(REALSXP, 2*nvar)));
    dtemp = REAL(SET_VECTOR_ELT(rlist, 1, allocMatrix(REALSXP, 2*nvar, 2*nvar)));
    imat = dmatrix(dtemp, 2*nvar, 2*nvar);  /* information matrix */
    
    dtemp = REAL(SET_VECTOR_ELT(rlist, 2, allocMatrix(REALSXP, nevent, nvar)));
    schoen = dmatrix(dtemp, nevent, nvar);  /* schoenfeld residuals */

    used = imatrix(INTEGER(SET_VECTOR_ELT(rlist, 3, 
					  allocMatrix(INTSXP, nstrat, nvar))),
		   nstrat, nvar);
    
    /* scratch vectors */
    a = (double *) R_alloc(2*nvar*nvar + 2*nvar, sizeof(double));
    a2 = a + nvar;
    cmat = dmatrix(a2+ nvar,   nvar, nvar);
    cmat2= dmatrix(a2+ nvar +nvar*nvar, nvar, nvar);

    /* count the number of events, per covariate per strata 
    ** if a covariate is constant in the stratum it gets a 0, otherwise
    **  the number of events in the stratum.  The algorithm is fairly
    **  simple but the code is ugly.  Count the number of deaths, then at
    **  the end of each stratum set used[][] to the number of events or
    **  to zero for each covariate.  
    */
    k = 0;  /* first obs of the stratum */
    ndead=0;
    cstrat = strata[sort[0]];
    for (i=0; i<nused; i++) {
	person= sort[i];
	if (cstrat == strata[person]) ndead += status[person];
	else { /* end of a stratum */
	    for (j=0; j<nvar; j++) {
		used[j][cstrat] =0;   /* start pessimistic */
		for (kk =k; kk<i; kk++) {
		    person = sort[kk];
		    if (covar[j][person] != covar[j][sort[k]]) {
			used[j][cstrat] = ndead;
			break;
		    }
		}
	    }
	    ndead= status[sort[i]];
	    k = i;
	    cstrat = strata[sort[i]];
	}
    }
    /* Deal with the last strata */
    for (j=0; j<nvar; j++) {
	used[j][cstrat] =0;   /* start pessimistic */
	for (kk =k; kk<i; kk++) {
	    person = sort[kk];
	    if (covar[j][person] != covar[j][sort[k]]) {
		used[j][cstrat] = ndead;
		break;
	    }
	}
    }
 	
    /*
    **	Recenter the X matrix to make the variance computation more stable
    */
    for (i=0; i<nvar; i++) {
	tmean =0;
	for (j=0; j<nused; j++) {
	    tmean += covar[i][j];
	}
	tmean /= nused;
	for (j=0; j<nused; j++) covar[i][j] -= tmean;
    }	


    /* zero variables */
    for (i=0; i<2*nvar; i++) {
	u[i] =0;
	for (j=0; j<2*nvar; j++) imat[i][j] =0 ;
	}
    for (i=0; i<nvar; i++) {
	a[i] =0;
	a2[i] =0;
	for (j=0; j<nvar; j++) {
	    cmat[i][j] =0;
	    cmat2[i][j] =0;
	}
    }

    /*
    ** Compute first and second derivatives
    */
    cstrat = -1; /* will not match any data point */
    ip = nused -1;
    while(ip >=0) {  /* backwards in time */
	person = sort[ip];
	if (strata[person] != cstrat) {
	    cstrat= strata[person];
	    denom = 0;
	    for (i=0; i<nvar; i++) {
		a[i] = 0;
		for (j=0; j<nvar; j++) cmat[i][j] = 0;
		}
	    }

	dtime = xtime[person];
	timewt = gt[person];   /* time weight for this event time */
	ndead =0; /*number of deaths at this time point */
	deadwt =0;  /* sum of weights for the deaths */
	denom2=0;  /* sum of weighted risks for the deaths */
	for (; ip>=0; ip--) {
	    person = sort[ip];
	    if (xtime[person] != dtime || strata[person] !=cstrat) break;
	    /* walk through the this set of tied times */
	    risk = exp(eta[person]) * weights[person];
	    if (status[person] ==0) {
		denom += risk;
		/* a contains weighted sums of x, cmat sums of squares */
		for (i=0; i<nvar; i++) {
		    a[i] += risk*covar[i][person];
		    for (j=0; j<=i; j++)
			cmat[i][j] += risk*covar[i][person]*covar[j][person];
	        }
	    }	
	    else {
		ndead++;
		deadwt += weights[person];
		denom2 += risk;
		nevent--;
		for (i=0; i<nvar; i++) {
		    schoen[i][nevent] = covar[i][person];
		    u[i] += weights[person]*covar[i][person];
		    u[i+nvar] += timewt* weights[person]*covar[i][person];
		    a2[i] +=  risk*covar[i][person];
		    for (j=0; j<=i; j++)
			cmat2[i][j] += risk*covar[i][person]*covar[j][person];
		        }
	    }
	}

	if (ndead >0) {  /* we need to add to the main terms */
	    if (method==0) { /* Breslow */
		denom += denom2;
	   
		for (i=0; i<nvar; i++) {
		    a[i] += a2[i];
		    temp2= a[i]/ denom;  /* mean */
		    u[i] -=  deadwt* temp2;
		    u[i+nvar] -= timewt* deadwt* temp2;
		    for (j=0; j<=i; j++) {
			cmat[i][j] += cmat2[i][j];
			temp = deadwt*(cmat[i][j] - temp2*a[j])/denom;
			imat[j][i] += temp;
			imat[j][i+nvar] += temp* timewt;
			imat[j+nvar][i+nvar] += temp*timewt * timewt;
			}
		    for (j=0; j<ndead; j++)
			schoen[i][nevent+j] -= temp2;
		    }
		}
	    else { /* Efron */
		/*
		** If there are 3 deaths we have 3 terms: in the first the
		**  three deaths are all in, in the second they are 2/3
		**  in the sums, and in the last 1/3 in the sum.  Let k go
		**  1 to ndead: we sequentially add a2/ndead and cmat2/ndead
		**  and efron_wt/ndead to the totals.
		*/
		wtave = deadwt/ndead;

		/* compute the mean of the means, for Schoenfeld resids */
		for (i=0; i<nvar; i++) {
		    tmean =0;
		    for (k=0; k<ndead; k++) {
			temp = (k+1)/ndead;
			tmean += (a[i] + a2[i]*temp)/(denom + denom2*temp);
			}
		    for (j=0; j<ndead; j++)
			schoen[i][nevent+j] -= tmean/ndead;
		}

		/* now compute U and imat */
		for (k=0; k<ndead; k++) {
		    denom += denom2/ndead;
		    for (i=0; i<nvar; i++) {
			a[i] += a2[i]/ndead;
			temp2 = a[i]/denom;
			u[i] -= wtave *temp2;
			u[i+nvar] -= timewt* wtave *temp2;
			for (j=0; j<=i; j++) {
			    cmat[i][j] += cmat2[i][j]/ndead;
			    temp =  wtave*(cmat[i][j] - temp2*a[j])/denom;
			    imat[j][i] += temp;
			    imat[j][i+nvar] += timewt*temp;
			    imat[j+nvar][i+nvar] += timewt * timewt *temp;
			}	
		    }
		}
	    }	
	    for (i=0; i<nvar; i++) {
		a2[i]=0;
		for (j=0; j<nvar; j++) cmat2[i][j]=0;
	    }
	}
    } /* end  of accumulation loop */

    /* fill in the rest of the information matrix */
    for (i=0; i<nvar; i++) {
	for (j=0; j<i; j++) {
	    imat[i][j] = imat[j][i];
	    imat[i][j+nvar]= imat[j][i+nvar];
	    imat[i+nvar][j+nvar] = imat[j+nvar][i+nvar];
	}
    }	
    for (i=0; i<nvar; i++) {
	for (j=0; j<nvar; j++)  /* upper right */
	    imat[i+nvar][j] = imat[j][i+nvar];
    }	

    UNPROTECT(nprotect);
    return(rlist);
}	
