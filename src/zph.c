/* 
** a 1 step variant of coxph5 for the zph test.
** gt = a function of time, with one element per unique death time
**   (within strata)
**  only compute the information needed for the score test,
**  which is the vector \sum g(t)* (x_i - xbar) of p elements and the 
**  full 2p by 2p information matrix
*/

#include <math.h>
#include "survS.h"
#include "survproto.h"

SEXP zph1(SEXP gt2,    SEXP y2, 
	  SEXP covar2, SEXP eta2,  SEXP weights2,
	  SEXP strata2,SEXP method2, SEXP sort2) {

    int i,j, k, person, ip;
    int cstrat;   /* current stratum*/
    double temp, temp2;
    double *dtemp;
    double *a, *a2, **cmat, **cmat2;
    double nrisk, denom, dtime, ndead, deadwt, denom2;
    double wtave, risk;

    /* scalar input arguments */
    int     nused, nvar;
    int     method, ntime;

    /* input vectors */
    double *xtime, *gt, *eta, *weights, *status;
    int    *strata, *sort;
    double **covar, **imat;

    /* returned objects */
    SEXP rlist;
    double *u;
    static const char *rnames[]={"u", "imat", ""};

    /* get local copies of input args */
    nused = nrows(y2);
    nvar  = ncols(covar2);
    method = asInteger(method2);
    ntime = LENGTH(gt2);

    xtime =   REAL(y2);
    status =  xtime + nused;
    weights = REAL(weights2);
    strata = INTEGER(strata2);
    gt =     REAL(gt2);
    eta =    REAL(eta2);
    sort=    INTEGER(sort2);

    /*
    **  Set up the ragged arrays and scratch space
    */
    PROTECT(rlist = mkNamed(VECSXP, rnames));
    covar= dmatrix(REAL(covar2), nused, nvar);
    u = REAL(SET_VECTOR_ELT(rlist, 0, allocVector(REALSXP, 2*nvar)));
    dtemp = REAL(SET_VECTOR_ELT(rlist, 1, allocMatrix(REALSXP, 2*nvar, 2*nvar)));
    imat = dmatrix(dtemp, 2*nvar, 2*nvar);
    a = (double *) R_alloc(2*nvar*nvar + 2*nvar, sizeof(double));
    a2 = a + nvar;
    cmat = dmatrix(a2+ nvar,   nvar, nvar);
    cmat2= dmatrix(a2+ nvar +nvar*nvar, nvar, nvar);


    /*
    ** Compute first and second derivatives
    */
    cstrat = -1; /* will not match any data point */

    for (i=0; i<2*nvar; i++) {
	u[i] =0;
	for (j=0; j<2*nvar; j++) imat[i][j] =0 ;
	}
    for (i=0; i<nvar; i++) {
	a2[i] =0;
	for (j=0; j<nvar; j++) cmat2[i][j] =0;
	}

    for (ip=0; ip<nused; ip++) {
	person = sort[i];
	if (strata[person] != cstrat) {
	    cstrat= strata[person];
	    nrisk =0 ;  
	    denom = 0;
	    for (i=0; i<nvar; i++) {
		a[i] = 0;
		for (j=0; j<nvar; j++) cmat[i][j] = 0;
		}
	    }

	dtime = xtime[person];
	ndead =0; /*number of deaths at this time point */
	deadwt =0;  /* sum of weights for the deaths */
	denom2=0;  /* sum of weighted risks for the deaths */
	for(; ip<nused; ip++) {
	    person= sort[ip] -1;
	    if (xtime[person] != dtime || strata[person] != cstrat) break;

	    /* walk through the this set of tied times */
	    nrisk++;
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
		ntime--;
		if (ntime <0) error("ntime error in zph.c");
		ndead++;
		deadwt += weights[person];
		denom2 += risk;

		for (i=0; i<nvar; i++) {
		    u[i] += weights[person]*covar[i][person];
		    u[i+nvar] += gt[ntime]* weights[person]*covar[i][person];
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
		    u[i+nvar] -= gt[ntime]* deadwt* temp2;
		    for (j=0; j<=i; j++) {
			cmat[i][j] += cmat2[i][j];
			temp = deadwt*(cmat[i][j] - temp2*a[j])/denom;
			imat[j][i] += temp;
			imat[j][i+nvar] += temp* gt[ntime];
			imat[j+nvar][i+nvar] += temp*gt[ntime]*gt[ntime];
			}
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
		for (k=0; k<ndead; k++) {
		    denom += denom2/ndead;
		    for (i=0; i<nvar; i++) {
			a[i] += a2[i]/ndead;
			temp2 = a[i]/denom;
			u[i] -= wtave *temp2;
			u[i+nvar] -= gt[ntime]* wtave *temp2;
			for (j=0; j<=i; j++) {
			    cmat[i][j] += cmat2[i][j]/ndead;
			    temp =  wtave*(cmat[i][j] - temp2*a[j])/denom;
			    imat[j][i] += temp;
			    imat[j][i+nvar] += gt[ntime]*temp;
			    imat[j+nvar][i+nvar] += gt[ntime]*gt[ntime]*temp;
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
    UNPROTECT(1);
    return(rlist);
}	
