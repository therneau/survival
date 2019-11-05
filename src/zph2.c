/*
** The code for start-stop data
**  this is the essentially the first iteration of agfit4.c
*/
#include <math.h>
#include "survS.h"
#include "survproto.h"

SEXP zph2(SEXP gt2,    SEXP y2, 
	  SEXP covar2, SEXP eta2,  SEXP weights2,
	  SEXP strata2,SEXP method2, SEXP sort12, SEXP sort22) {

    int i,j, k, kk, person, p, p1;
    int cstrat;   /* current stratum*/
    int indx1, nrisk;
    double temp, temp2, tmean, etasum;
    double *dtemp, timewt=1;
    double *a, *a2, **cmat, **cmat2;
    int *keep;
    double denom, dtime=0, ndead, denom2;
    double risk, meanwt;

    /* scalar input arguments and counts*/
    int  nused, nvar, nevent, nstrat;
    int  method;

    /* input vectors */
    double *start, *tstop, *gt, *eta, *weights, *status;
    int    *strata, *sort1, *sort2;
    double **covar;

    /* returned objects */
    SEXP rlist;
    double *u, **imat, **schoen;
    int **used;
    static const char *rnames[]={"u", "imat", "schoen", "used", ""};

    /* get local copies of input args */
    nused = nrows(y2);
    nvar  = ncols(covar2);
    method = asInteger(method2);

    start =   REAL(y2);
    tstop =   start + nused;
    status =  tstop + nused;
    weights = REAL(weights2);
    strata = INTEGER(strata2);
    gt =     REAL(gt2);
    eta =    REAL(eta2);
    sort1=   INTEGER(sort12);
    sort2=   INTEGER(sort22);
 
   /* 
    ** count up the number of events and the number of strata
    */
    nevent =0;
    nstrat =1;
    for (i=0; i< nused; i++) {
	nevent += status[i];
	if (strata[i] >= nstrat) nstrat= strata[i] + 1;
    }

    /*
    **  Set up the ragged arrays and scratch space
    */
    if (MAYBE_REFERENCED(covar2))
	covar = dmatrix(REAL(duplicate(covar2)), nused, nvar);
    else covar = dmatrix(REAL(covar2), nused, nvar);

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
    keep = (int *) R_alloc(nused, sizeof(int));

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
	for (j=0; j<nstrat; j++) used[i][j] =0;
	}
    for (i=0; i<nused; i++) keep[i] =0;

    /* count the number of events, per covariate per strata 
    ** if a covariate is constant in the stratum it gets a 0, otherwise
    **  the number of events in the stratum.  The algorithm is fairly
    **  simple but the code is ugly.  Count the number of deaths, then at
    **  the end of each stratum set used[][] to the number of events or
    **  to zero for each covariate.  
    */
    k = 0;  /* first obs of the stratum */
    ndead=0;
    cstrat = strata[sort2[0]];
    for (i=0; i<nused; i++) {
	person= sort2[i];
	if (cstrat == strata[person]) ndead += status[person];
	else { /* end of a stratum */
	    p1 = sort2[k];   /* first obs of the stratum */
	    for (j=0; j<nvar; j++) {
		used[j][cstrat] =0;   /* start pessimistic */
		for (kk =k; kk<i; kk++) {
		    person = sort2[kk];
		    if (covar[j][person] != covar[j][p1]) {
			used[j][cstrat] = ndead;
			break;
		    }
		}
	    }
	    k = i;
	    ndead = status[sort2[i]];
	    cstrat = strata[sort2[i]];
	}
    }

    /* Deal with the last strata */
    p1 = sort2[k];
    for (j=0; j<nvar; j++) {
	used[j][cstrat] =0;   /* start pessimistic */
	for (kk =k; kk<i; kk++) {
	    person = sort2[kk];
	    if (covar[j][person] != covar[j][p1]) {
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
	for (j=0; j<nused; j++) tmean += covar[i][j];
	tmean /= nused;
	for (j=0; j<nused; j++) covar[i][j] -= tmean;
    }	

    /*
    **  'person' walks through the the data from 1 to n,
    **     sort1[0] points to the largest stop time, sort1[1] the next, ...
    **  'dtime' is a scratch variable holding the time of current interest
    **  'indx1' walks through the start times.  
    */
    person =0;
    indx1 =0;
    denom=0;
    nrisk=0;
    etasum =0;
    cstrat = -1;
    while (person < nused) {
        /* find the next death time */
        for (k=person; k< nused; k++) {
            if (strata[sort2[k]] != cstrat) {
                /* hit a new stratum; reset temporary sums */
                cstrat = strata[sort2[k]];
                denom = 0;
                nrisk = 0;
                etasum =0;
                for (i=0; i<nvar; i++) {
                    a[i] =0;
                    for (j=0; j<nvar; j++) cmat[i][j] =0;
                }
                person =k;  /* skip to end of stratum */
                indx1  =k; 
            }
            p = sort2[k];
            if (status[p] == 1) {
                dtime = tstop[p];
		timewt = gt[p];   /* time weight for this event time */
                break;
            }
        }
        if (k == nused) person =k;  /* no more ndead to be processed */
        else {
            /* remove any subjects no longer at risk */
            /*
            ** subtract out the subjects whose start time is to the right
            ** If everyone is removed reset the totals to zero. 
            */
            for (; indx1<nused; indx1++) {
                p1 = sort1[indx1];
                if (strata[p1] != cstrat || start[p1] < dtime) break;
                if (keep[p1] == 0) continue;  /* skip any never-at-risk rows */
                nrisk--;
                if (nrisk ==0) {
                    etasum =0;
                    denom =0;
                    for (i=0; i<nvar; i++) {
                        a[i] =0;
                        for (j=0; j<=i; j++) cmat[i][j] =0;
                    }
                }
                else {
                    etasum -= eta[p1];
                    risk = exp(eta[p1]) * weights[p1];
                    denom -= risk;
                    for (i=0; i<nvar; i++) {
                        a[i] -= risk*covar[i][p1];
                        for (j=0; j<=i; j++)
                            cmat[i][j] -= risk*covar[i][p1]*covar[j][p1];
                    }
                }
                if (fabs(etasum/nrisk) > 200) {  
                    temp = etasum/nrisk;
                    for (i=0; i<nused; i++) eta[i] -= temp;
                    temp = exp(-temp);
                    denom *= temp;
                    for (i=0; i<nvar; i++) {
                        a[i] *= temp;
                        for (j=0; j<nvar; j++) {
                            cmat[i][j]*= temp;
                        }
                    }
                    etasum =0;
                }
            }

            /* 
            ** add any new subjects who are at risk 
            ** denom2, a2, cmat2, meanwt and ndead count only the ndead
            */
            denom2= 0;
            ndead=0;  meanwt=0;  
            for (i=0; i<nvar; i++) {
                a2[i]=0;
                for (j=0; j<nvar; j++) {
                    cmat2[i][j]=0;
                }
            }
            
            for (; person<nused; person++) {
                p = sort2[person];
                if (strata[p] != cstrat || tstop[p] < dtime) break; /* no more to add */
                risk = exp(eta[p]) * weights[p];

                if (status[p] ==1) {
		    nevent--;
                    keep[p] =1;
                    nrisk++;
                    etasum += eta[p];
                    ndead++;
                    denom2 += risk;
                    meanwt += weights[p];
                    for (i=0; i<nvar; i++) {
                        u[i] += weights[p] * covar[i][p];
			u[i+nvar] += weights[p]*covar[i][p] * timewt;
                        a2[i]+= risk*covar[i][p];
			schoen[i][nevent] = covar[i][p];
                        for (j=0; j<=i; j++)
                            cmat2[i][j] += risk*covar[i][p]*covar[j][p];
                    }
                }
                else if (start[p] < dtime) {
                    keep[p] =1;
                    nrisk++;
                    etasum += eta[p];
                    denom += risk;
                    for (i=0; i<nvar; i++) {
                        a[i] += risk*covar[i][p];
                        for (j=0; j<=i; j++)
                            cmat[i][j] += risk*covar[i][p]*covar[j][p];
                    }
                } 
            }
            /*
            ** Add results into u and imat for all events at this time point
            */
            if (method==0 || ndead ==1) { /*Breslow */
                denom += denom2;
                for (i=0; i<nvar; i++) {
                    a[i] += a2[i];
                    temp = a[i]/denom;   /*mean covariate at this time */
                    u[i] -= meanwt*temp;
		    u[i+nvar] -= meanwt*temp * timewt;
                    for (j=0; j<=i; j++) {
                        cmat[i][j] += cmat2[i][j];
			temp2 = meanwt*((cmat[i][j]- temp*a[j])/denom);
                        imat[j][i] += temp2;
			imat[j][i+nvar] += temp2* timewt;
			imat[j+nvar][i+nvar] += temp2*timewt * timewt;
                    }
		    for (j=0; j<ndead; j++)
			schoen[i][nevent+j] -= temp;
                }
            }
            else {
                meanwt /= ndead;
		/* compute the mean x, for Schoenfeld residuals */
		for (i=0; i<nvar; i++) {
		    tmean =0;
		    for (k=0; k<ndead; k++) {
			temp = (k+1)/ndead;
			tmean += (a[i] + a2[i]*temp)/(denom + denom2*temp);
			}
		    for (j=0; j<ndead; j++)
			schoen[i][nevent+j] -= tmean/ndead;
		}

                for (k=0; k<ndead; k++) {
                    denom += denom2/ndead;
		    for (i=0; i<nvar; i++) {
                        a[i] += a2[i]/ndead;
                        temp = a[i]/denom;
                        u[i] -= meanwt*temp;
			u[i+nvar] -= meanwt*temp* timewt;
			for (j=0; j<=i; j++) {
                            cmat[i][j] += cmat2[i][j]/ndead;
			    temp2 = meanwt*((cmat[i][j]- temp*a[j])/denom);
                            imat[j][i] += temp2;
			    imat[j][i+nvar] += temp2* timewt;
			    imat[j+nvar][i+nvar] += temp2*timewt * timewt;
			} 
		    }
		}	
            }
            if (fabs(etasum/nrisk) > 200) {  
                temp = etasum/nrisk;
                for (i=0; i<nused; i++) eta[i] -= temp;
                temp = exp(-temp);
                denom *= temp;
                for (i=0; i<nvar; i++) {
                    a[i] *= temp;
                    for (j=0; j<nvar; j++) {
                        cmat[i][j]*= temp;
                    }
                }
                etasum =0;
            }
        }
    }   /* end  of accumulation loop */
    
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

    UNPROTECT(1);
    return(rlist);
}
