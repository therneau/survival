/* Automatically generated from the noweb directory */
#include <math.h>
#include "survS.h" 
#include "survproto.h"

SEXP agfit4(SEXP nused2, SEXP surv2,      SEXP covar2,    SEXP strata2,
            SEXP weights2,   SEXP offset2,   SEXP ibeta2,
            SEXP sort12,     SEXP sort22,    SEXP method2,
            SEXP maxiter2,   SEXP  eps2,     SEXP tolerance2,
            SEXP doscale2) { 
                
    int i,j,k, person;
    int indx1, istrat, p, p1;
    int nrisk, nr;
    int nused, nvar;
    int rank=0, rank2, fail;  /* =0 to keep -Wall happy */
   
    double **covar, **cmat, **imat;  /*ragged array versions*/
    double *a, *oldbeta;
    double *scale;
    double *a2, **cmat2;
    double *eta;
    double  denom, zbeta, risk;
    double  dtime =0;  /* initial value to stop a -Wall message */
    double  temp, temp2;
    double  newlk =0;
    int  halving;    /*are we doing step halving at the moment? */
    double  tol_chol, eps;
    double  meanwt;
    int deaths;
    double denom2, etasum;
 
    /* inputs */
    double *start, *tstop, *event;
    double *weights, *offset;
    int *sort1, *sort2, maxiter;
    int *strata;
    double method;  /* saving this as double forces some double arithmetic */
    int doscale;

    /* returned objects */
    SEXP imat2, beta2, u2, loglik2;
    double *beta, *u, *loglik;
    SEXP sctest2, flag2, iter2;
    double *sctest;
    int *flag, *iter;
    SEXP rlist;
    static const char *outnames[]={"coef", "u", "imat", "loglik",
                                   "sctest", "flag", "iter", ""};
    int nprotect;  /* number of protect calls I have issued */

    /* get sizes and constants */
    nused = asInteger(nused2);
    nvar  = ncols(covar2);
    nr    = nrows(covar2);  /*nr = number of rows, nused = how many we use */
    method= asInteger(method2);
    eps   = asReal(eps2);
    tol_chol = asReal(tolerance2);
    maxiter = asInteger(maxiter2);
    doscale = asInteger(doscale2);
  
    /* input arguments */
    start = REAL(surv2);
    tstop  = start + nr;
    event = tstop + nr;
    weights = REAL(weights2);
    offset = REAL(offset2);
    sort1  = INTEGER(sort12);
    sort2  = INTEGER(sort22);
    strata = INTEGER(strata2);

    /*
    ** scratch space
    **  nvar: a, a2, oldbeta, scale
    **  nvar*nvar: cmat, cmat2
    **  nr:  eta
    */
    eta = (double *) R_alloc(nr + 4*nvar + 2*nvar*nvar, sizeof(double));
    a = eta + nr;
    a2= a + nvar;
    scale  = a2 + nvar;
    oldbeta = scale + nvar;
            
    /*
    **  Set up the ragged arrays
    **  covar2 might not need to be duplicated, even though
    **  we are going to modify it, due to the way this routine was
    **  was called.  But check
    */
    PROTECT(imat2 = allocMatrix(REALSXP, nvar, nvar));
    nprotect =1;
    if (MAYBE_REFERENCED(covar2)) {
        PROTECT(covar2 = duplicate(covar2)); 
        nprotect++;
        }
    covar= dmatrix(REAL(covar2), nr, nvar);
    imat = dmatrix(REAL(imat2),  nvar, nvar);
    cmat = dmatrix(oldbeta+ nvar,   nvar, nvar);
    cmat2= dmatrix(oldbeta+ nvar + nvar*nvar, nvar, nvar);

    /*
    ** create the output structures
    */
    PROTECT(rlist = mkNamed(VECSXP, outnames));
    nprotect++;
    beta2 = SET_VECTOR_ELT(rlist, 0, duplicate(ibeta2));
    beta  = REAL(beta2);
    u2 =    SET_VECTOR_ELT(rlist, 1, allocVector(REALSXP, nvar));
    u = REAL(u2);

    SET_VECTOR_ELT(rlist, 2, imat2);
    loglik2 = SET_VECTOR_ELT(rlist, 3, allocVector(REALSXP, 2)); 
    loglik  = REAL(loglik2);

    sctest2 = SET_VECTOR_ELT(rlist, 4, allocVector(REALSXP, 1));
    sctest =  REAL(sctest2);
    flag2  =  SET_VECTOR_ELT(rlist, 5, allocVector(INTSXP, 4));
    flag   =  INTEGER(flag2);
    for (i=0; i<4; i++) flag[i]=0;

    iter2  =  SET_VECTOR_ELT(rlist, 6, allocVector(INTSXP, 1));
    iter = INTEGER(iter2);
                
    /*
    ** Subtract the mean from each covar, as this makes the variance
    **  computation much more stable.  The mean is taken per stratum,
    **  the scaling is overall.
    */
    for (i=0; i<nvar; i++) {
        istrat = strata[sort2[0]];  /* the current stratum */
        k = 0;                      /* first obs of current one */
        temp =0;  temp2=0;
        for (person=0; person< nused; person++) {
            p = sort2[person];
            if (strata[p] == istrat) {
                temp += weights[p] * covar[i][p];
                temp2 += weights[p];
            }
            else {  /* new stratum */
                temp /= temp2;  /* mean for this covariate, this strata */
                for (; k< person; k++) covar[i][sort2[k]] -=temp;
                temp =0;  temp2=0;
                istrat = strata[p];
            }
        temp /= temp2;  /* mean for last stratum */
        for (; k< nused; k++) covar[i][sort2[k]] -= temp;
        }
        if (doscale ==1) { /* also scale the regression */
            /* this cannot be done per stratum */
            temp =0;
            temp2 =0;
            for (person=0; person<nused; person++) {
                p = sort2[person];
                temp += weights[p] * fabs(covar[i][p]);
                temp2 += weights[p];
                }
            if (temp >0) temp = temp2/temp;  /* 1/scale */
            else temp = 1.0;  /* rare case of a constant covariate */
            scale[i] = temp;
            for (person=0; person<nused; person++) {
                covar[i][sort2[person]] *= temp;
            }
        }
    }
 
    if (doscale ==1) {
        for (i=0; i<nvar; i++) beta[i] /= scale[i]; /* rescale initial betas */
        }
    else {for (i=0; i<nvar; i++) scale[i] = 1.0;}
             
    /* main loop */
    halving =0 ;             /* =1 when in the midst of "step halving" */
    fail =0;
    for (*iter=0; *iter<= maxiter; (*iter)++) {
        R_CheckUserInterrupt();  /* be polite -- did the user hit cntrl-C? */
        for (person=0; person<nused; person++) {
            p = sort2[person];
            zbeta = 0;      /* form the term beta*z   (vector mult) */
            for (i=0; i<nvar; i++)
                zbeta += beta[i]*covar[i][p];
            eta[p] = zbeta + offset[p];
        }

        /*
        **  'person' walks through the the data from 1 to nused,
        **     sort1[0] points to the largest stop time, sort1[1] the next, ...
        **  'dtime' is a scratch variable holding the time of current interest
        **  'indx1' walks through the start times.  
        */
        newlk =0;
        for (i=0; i<nvar; i++) {
            u[i] =0;
            for (j=0; j<nvar; j++) imat[i][j] =0;
        }
        person =0;
        indx1 =0;

        /* this next set is rezeroed at the start of each stratum */
        denom=0;
        nrisk=0;
        etasum =0;
        for (i=0; i<nvar; i++) {
            a[i] =0;
            for (j=0; j<nvar; j++) cmat[i][j] =0;
        }
        /* end of the per-stratum set */

        istrat = strata[sort2[0]];  /* initial stratum */
        while (person < nused) {
            /* find the next death time */
            for (k=person; k< nused; k++) {
                p = sort2[k];
                if (strata[p] != istrat) {
                    /* hit a new stratum; reset temporary sums */
                    istrat= strata[p];
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

                if (event[p] == 1) {
                    dtime = tstop[p];
                    break;
                }
            }
            if (k == nused) break;  /* no more deaths to be processed */

            /* remove any subjects no longer at risk */
            /*
            ** subtract out the subjects whose start time is to the right
            ** If everyone is removed reset the totals to zero.  (This happens when
            ** the survSplit function is used, so it is worth checking).
            */
            for (; indx1<nused; indx1++) {
                p1 = sort1[indx1];
                if (start[p1] < dtime || strata[p1] != istrat) break;
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
                /* 
                ** We must avoid overflow in the exp function (~750 on Intel)
                ** and want to act well before that, but not take action very often.  
                ** One of the case-cohort papers suggests an offset of -100 meaning
                ** that etas of 50-100 can occur in "ok" data, so make it larger 
                ** than this.
                ** If the range of eta is more then log(1e16) = 37 then the data is
                **  hopeless: some observations will have effectively 0 weight.  Keeping
                **  the mean sensible suffices to keep the max in check for all other
                *   data sets.
                */
                if (fabs(etasum/nrisk) > 200) {  
                    flag[1]++;  /* a count, for debugging/profiling purposes */
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
            ** denom2, a2, cmat2, meanwt and deaths count only the deaths
            */
            denom2= 0;
            meanwt =0;
            deaths=0;    
            for (i=0; i<nvar; i++) {
                a2[i]=0;
                for (j=0; j<nvar; j++) {
                    cmat2[i][j]=0;
                }
            }
            
            for (; person <nused; person++) {
                p = sort2[person];
                if (strata[p] != istrat || tstop[p] < dtime) break;/*no more to add*/
                risk = exp(eta[p]) * weights[p];
                
                if (event[p] ==1 ){
                    nrisk++;
                    etasum += eta[p];
                    deaths++;
                    denom2 += risk;
                    meanwt += weights[p];
                    newlk += weights[p]* eta[p];
                    for (i=0; i<nvar; i++) {
                        u[i] += weights[p] * covar[i][p];
                        a2[i]+= risk*covar[i][p];
                        for (j=0; j<=i; j++)
                            cmat2[i][j] += risk*covar[i][p]*covar[j][p];
                    }
                }
                else {
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
            if (method==0 || deaths ==1) { /*Breslow */
                denom += denom2;
                newlk -= meanwt*log(denom);  /* sum of death weights*/ 
                for (i=0; i<nvar; i++) {
                    a[i] += a2[i];
                    temp = a[i]/denom;   /*mean covariate at this time */
                    u[i] -= meanwt*temp;
                    for (j=0; j<=i; j++) {
                        cmat[i][j] += cmat2[i][j];
                        imat[j][i] += meanwt*((cmat[i][j]- temp*a[j])/denom);
                    }
                }
            }
            else {
                meanwt /= deaths;
                for (k=0; k<deaths; k++) {
                    denom += denom2/deaths;
                    newlk -= meanwt*log(denom);
                    for (i=0; i<nvar; i++) {
                        a[i] += a2[i]/deaths;
                        temp = a[i]/denom;
                        u[i] -= meanwt*temp;
                        for (j=0; j<=i; j++) {
                            cmat[i][j] += cmat2[i][j]/deaths;
                            imat[j][i] += meanwt*((cmat[i][j]- temp*a[j])/denom);
                        }
                        }
                }
            }
            /* 
            ** We must avoid overflow in the exp function (~750 on Intel)
            ** and want to act well before that, but not take action very often.  
            ** One of the case-cohort papers suggests an offset of -100 meaning
            ** that etas of 50-100 can occur in "ok" data, so make it larger 
            ** than this.
            ** If the range of eta is more then log(1e16) = 37 then the data is
            **  hopeless: some observations will have effectively 0 weight.  Keeping
            **  the mean sensible suffices to keep the max in check for all other
            *   data sets.
            */
            if (fabs(etasum/nrisk) > 200) {  
                flag[1]++;  /* a count, for debugging/profiling purposes */
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
        }   /* end  of accumulation loop */

        if (*iter==0) {
            loglik[0] = newlk;
            loglik[1] = newlk;
            /* compute the score test, but don't corrupt u */
            for (i=0; i<nvar; i++) a[i] = u[i];
            rank = cholesky2(imat, nvar, tol_chol);
            chsolve2(imat,nvar,a);        /* a replaced by  u *inverse(i) */
            *sctest=0;
            for (i=0; i<nvar; i++) {
               *sctest +=  u[i]*a[i];
            }
            if (maxiter==0) break;
            for (i=0; i<nvar; i++) {
                  oldbeta[i] = beta[i];
                beta[i] += a[i];
            }        
        }
        else { 
            if (*iter ==1) {
                fail = isnan(newlk) + isinf(newlk);
                /* it almost takes malice to give a starting estimate with infinite
                **  loglik.  But if so, just give up now */
                 if (fail>0) break;
            }

            fail =0;
            for (i=0; i<nvar; i++) 
                if (isfinite(imat[i][i]) ==0) fail++;
            rank2 = cholesky2(imat, nvar, tol_chol);
            fail = fail + isnan(newlk) + isinf(newlk) + abs(rank-rank2);
     
            if (fail ==0 && halving ==0 &&
                fabs(1-(loglik[1]/newlk)) <= eps) break;  /* success! */

            if (*iter == maxiter) { /* failed to converge */
               flag[3] = 1;  
               if (maxiter>1 && ((newlk -loglik[1])/ loglik[1]) < -eps) {
                   /* 
                   ** "Once more unto the breach, dear friends, once more; ..."
                   **The last iteration above was worse than one of the earlier ones,
                   **  by more than roundoff error.  
                   ** We need to use beta and imat at the last good value, not the
                   **  last attempted value. We have tossed the old imat away, so 
                   **  recompute it.
                   ** It will happen very rarely that we run out of iterations, and
                   **  even less often that the Newton-Raphson is getting worse.
                   */
                   for (i=0; i<nvar; i++) beta[i] = oldbeta[i];
                   for (person=0; person<nused; person++) {
                       p = sort2[person];
                       zbeta = 0;      /* form the term beta*z   (vector mult) */
                       for (i=0; i<nvar; i++)
                           zbeta += beta[i]*covar[i][p];
                       eta[p] = zbeta + offset[p];
                   }

                   /*
                   **  'person' walks through the the data from 1 to nused,
                   **     sort1[0] points to the largest stop time, sort1[1] the next, ...
                   **  'dtime' is a scratch variable holding the time of current interest
                   **  'indx1' walks through the start times.  
                   */
                   newlk =0;
                   for (i=0; i<nvar; i++) {
                       u[i] =0;
                       for (j=0; j<nvar; j++) imat[i][j] =0;
                   }
                   person =0;
                   indx1 =0;

                   /* this next set is rezeroed at the start of each stratum */
                   denom=0;
                   nrisk=0;
                   etasum =0;
                   for (i=0; i<nvar; i++) {
                       a[i] =0;
                       for (j=0; j<nvar; j++) cmat[i][j] =0;
                   }
                   /* end of the per-stratum set */

                   istrat = strata[sort2[0]];  /* initial stratum */
                   while (person < nused) {
                       /* find the next death time */
                       for (k=person; k< nused; k++) {
                           p = sort2[k];
                           if (strata[p] != istrat) {
                               /* hit a new stratum; reset temporary sums */
                               istrat= strata[p];
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

                           if (event[p] == 1) {
                               dtime = tstop[p];
                               break;
                           }
                       }
                       if (k == nused) break;  /* no more deaths to be processed */

                       /* remove any subjects no longer at risk */
                       /*
                       ** subtract out the subjects whose start time is to the right
                       ** If everyone is removed reset the totals to zero.  (This happens when
                       ** the survSplit function is used, so it is worth checking).
                       */
                       for (; indx1<nused; indx1++) {
                           p1 = sort1[indx1];
                           if (start[p1] < dtime || strata[p1] != istrat) break;
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
                           /* 
                           ** We must avoid overflow in the exp function (~750 on Intel)
                           ** and want to act well before that, but not take action very often.  
                           ** One of the case-cohort papers suggests an offset of -100 meaning
                           ** that etas of 50-100 can occur in "ok" data, so make it larger 
                           ** than this.
                           ** If the range of eta is more then log(1e16) = 37 then the data is
                           **  hopeless: some observations will have effectively 0 weight.  Keeping
                           **  the mean sensible suffices to keep the max in check for all other
                           *   data sets.
                           */
                           if (fabs(etasum/nrisk) > 200) {  
                               flag[1]++;  /* a count, for debugging/profiling purposes */
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
                       ** denom2, a2, cmat2, meanwt and deaths count only the deaths
                       */
                       denom2= 0;
                       meanwt =0;
                       deaths=0;    
                       for (i=0; i<nvar; i++) {
                           a2[i]=0;
                           for (j=0; j<nvar; j++) {
                               cmat2[i][j]=0;
                           }
                       }
                       
                       for (; person <nused; person++) {
                           p = sort2[person];
                           if (strata[p] != istrat || tstop[p] < dtime) break;/*no more to add*/
                           risk = exp(eta[p]) * weights[p];
                           
                           if (event[p] ==1 ){
                               nrisk++;
                               etasum += eta[p];
                               deaths++;
                               denom2 += risk;
                               meanwt += weights[p];
                               newlk += weights[p]* eta[p];
                               for (i=0; i<nvar; i++) {
                                   u[i] += weights[p] * covar[i][p];
                                   a2[i]+= risk*covar[i][p];
                                   for (j=0; j<=i; j++)
                                       cmat2[i][j] += risk*covar[i][p]*covar[j][p];
                               }
                           }
                           else {
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
                       if (method==0 || deaths ==1) { /*Breslow */
                           denom += denom2;
                           newlk -= meanwt*log(denom);  /* sum of death weights*/ 
                           for (i=0; i<nvar; i++) {
                               a[i] += a2[i];
                               temp = a[i]/denom;   /*mean covariate at this time */
                               u[i] -= meanwt*temp;
                               for (j=0; j<=i; j++) {
                                   cmat[i][j] += cmat2[i][j];
                                   imat[j][i] += meanwt*((cmat[i][j]- temp*a[j])/denom);
                               }
                           }
                       }
                       else {
                           meanwt /= deaths;
                           for (k=0; k<deaths; k++) {
                               denom += denom2/deaths;
                               newlk -= meanwt*log(denom);
                               for (i=0; i<nvar; i++) {
                                   a[i] += a2[i]/deaths;
                                   temp = a[i]/denom;
                                   u[i] -= meanwt*temp;
                                   for (j=0; j<=i; j++) {
                                       cmat[i][j] += cmat2[i][j]/deaths;
                                       imat[j][i] += meanwt*((cmat[i][j]- temp*a[j])/denom);
                                   }
                                   }
                           }
                       }
                       /* 
                       ** We must avoid overflow in the exp function (~750 on Intel)
                       ** and want to act well before that, but not take action very often.  
                       ** One of the case-cohort papers suggests an offset of -100 meaning
                       ** that etas of 50-100 can occur in "ok" data, so make it larger 
                       ** than this.
                       ** If the range of eta is more then log(1e16) = 37 then the data is
                       **  hopeless: some observations will have effectively 0 weight.  Keeping
                       **  the mean sensible suffices to keep the max in check for all other
                       *   data sets.
                       */
                       if (fabs(etasum/nrisk) > 200) {  
                           flag[1]++;  /* a count, for debugging/profiling purposes */
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
                   }   /* end  of accumulation loop */
                   rank2 = cholesky2(imat, nvar, tol_chol);
                   }
               break;
            }
            
            if (fail >0 || newlk < loglik[1]) {
                /* 
                ** The routine has not made progress past the last good value.
                */
                halving++; flag[2]++;
                for (i=0; i<nvar; i++)
                    beta[i] = (oldbeta[i]*halving + beta[i]) /(halving +1.0);
            }
            else { 
                halving=0;
                loglik[1] = newlk;   /* best so far */  
                chsolve2(imat,nvar,u);
                for (i=0; i<nvar; i++) {
                    oldbeta[i] = beta[i];
                    beta[i] = beta[i] +  u[i];
                }
            }
        }
    } /*return for another iteration */

    flag[0] = rank; 
    loglik[1] = newlk;
    chinv2(imat, nvar);
    for (i=0; i<nvar; i++) {
        beta[i] *= scale[i];  /* return to original scale */
        u[i] /= scale[i];
        imat[i][i] *= scale[i] * scale[i];
        for (j=0; j<i; j++) {
            imat[j][i] *= scale[i] * scale[j];
            imat[i][j] = imat[j][i];
        }
    }
    UNPROTECT(nprotect);
    return(rlist);
}
