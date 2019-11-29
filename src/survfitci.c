/* Automatically generated from the noweb directory */
#include "survS.h"
#include "survproto.h"
#include <math.h>

SEXP survfitci(SEXP ftime2,  SEXP sort12,  SEXP sort22, SEXP ntime2,
                    SEXP status2, SEXP cstate2, SEXP wt2,  SEXP id2,
                    SEXP p2,      SEXP i02,     SEXP sefit2,
                    SEXP influence2,  SEXP inftime2,  SEXP time0) {   
    int i, j, k, kk;   /* generic loop indices */
    int ck, itime, eptr; /*specific indices */
    int itimej;
    double ctime;      /*current time of interest, in the main loop */
    int oldstate, newstate; /*when changing state */
    double lagtime;

    double temp, *temp2;  /* scratch double, and vector of length nstate */
    double *dptr;      /* reused in multiple contexts */
    double *p;         /* current prevalence vector */
    double **hmat;      /* hazard matrix at this time point */
    double **umat=0;     /* per subject leverage at this time point */
    int *atrisk;       /* 1 if the subject is currently at risk */
    int   *ns;         /* number curently in each state */
    int   *nev;        /* number of events at this time, by state */
    double *ws;        /* weighted count of number state */
    double *wtp;       /* case weights indexed by subject */
    double wevent;     /* weighted number of events at current time */
    int nstate;        /* number of states */
    int n, nperson;    /*number of obs, subjects*/
    double **chaz;     /* cumulative hazard matrix */

    /* pointers to the R variables */
    int *sort1, *sort2;  /*sort index for entry time, event time */
    double *entry,* etime;  /*entry time, event time */
    int ntime;          /* number of unique event time values */
    int *status;        /*0=censored, 1,2,... new states */
    int *cstate;        /* current state for each subject */
    int *dstate;        /* the next state, =cstate if not an event time */
    double *wt;         /* weight for each observation */
    double *i0;         /* initial influence */
    int *id;            /* for each obs, which subject is it */
    int sefit;
    int influence;
    double *inftime;
    int ninftime;
        
    /* returned objects */
    SEXP rlist;         /* the returned list and variable names of same */  
    const char *rnames[]= {"nrisk","nevent","ncensor", "p", 
                           "cumhaz", "std", "influence.pstate",
                           "influence.cumhaz", "influence.rmst", ""};
    SEXP setemp;
    double **pmat, **vmat=0, *cumhaz; /* =0 to silence -Wall warning */
    double *imat1=0;      /* returned leverage matrix, survival */
    double *imat3;        /* returned leverage for the RMST */
    int  *ncensor, **nrisk, **nevent;
    ntime= asInteger(ntime2);
    nperson = LENGTH(cstate2); /* number of unique subjects */
    n   = LENGTH(sort12);    /* number of observations in the data */
    PROTECT(cstate2 = duplicate(cstate2));
    cstate  = INTEGER(cstate2);
    entry= REAL(ftime2);
    etime= entry + n;
    sort1= INTEGER(sort12);
    sort2= INTEGER(sort22);
    status= INTEGER(status2);
    wt = REAL(wt2);
    id = INTEGER(id2);
    PROTECT(p2 = duplicate(p2));  /*copy of initial prevalence */
    p = REAL(p2);
    nstate = LENGTH(p2);  /* number of states */
    i0 = REAL(i02);
    sefit = asInteger(sefit2);
    influence = asInteger(influence2);
    inftime = REAL(inftime2);
    ninftime = LENGTH(inftime2);
    lagtime  = asReal(time0);

    /* allocate space for the output objects
    ** Ones that are put into a list do not need to be protected
    */
    PROTECT(rlist=mkNamed(VECSXP, rnames));
    setemp = SET_VECTOR_ELT(rlist, 0, allocMatrix(INTSXP, ntime, nstate));
    nrisk =  imatrix(INTEGER(setemp), ntime, nstate);  /* time by state */
    setemp = SET_VECTOR_ELT(rlist, 1, allocMatrix(INTSXP, ntime, nstate));
    nevent = imatrix(INTEGER(setemp), ntime, nstate);  /* time by state */
    setemp = SET_VECTOR_ELT(rlist, 2, allocVector(INTSXP, ntime));
    ncensor = INTEGER(setemp);  /* total at each time */
    setemp  = SET_VECTOR_ELT(rlist, 3, allocMatrix(REALSXP, ntime, nstate));
    pmat =   dmatrix(REAL(setemp), ntime, nstate);
    setemp = SET_VECTOR_ELT(rlist, 4, allocMatrix(REALSXP, nstate*nstate, ntime));
    cumhaz = REAL(setemp);

    if (sefit ==1) {
        setemp = SET_VECTOR_ELT(rlist, 5,  allocMatrix(REALSXP, ntime, nstate));
        vmat= dmatrix(REAL(setemp), ntime, nstate);
    }
    if (influence &0x1) {
        /* In R the max space is larger for a matrix than a vector.
        **  Use of allocVector(REALSXP, nperson * nstate * ninftime) may well 
        **    overflow the (integer) argument.  
        **  But alloc3DArray(REALSXP, nperson, nstate, ntime) is okay as long as
        **    each dimension is an integer.  There might not be enough memory to
        **    satisfy the request, but the request itself will be valid.
        **  Once allocated, I can treat setemp like a vector since imat is a pointer
        **   to double, a pointer is bigger than integer and will not overflow.
        */
        setemp = SET_VECTOR_ELT(rlist, 6, alloc3DArray(REALSXP, nperson, nstate, 
                                                        ninftime +1));
        imat1 = REAL(setemp);
    }
    if (influence &0x4) {
        setemp = SET_VECTOR_ELT(rlist, 8, alloc3DArray(REALSXP, nperson, nstate,
                                                        ninftime));
        imat3 = REAL(setemp);
    }  

    /* allocate space for scratch vectors */
    ws = (double *) R_alloc(2*nstate, sizeof(double)); /*weighted number in state */
    temp2 = ws + nstate;
    ns    = (int *) R_alloc(2*nstate, sizeof(int));
    nev   = ns + nstate;
    atrisk = (int *) R_alloc(2*nperson, sizeof(int));
    dstate = atrisk + nperson;
    wtp = (double *) R_alloc(nperson, sizeof(double));
    hmat = (double**) dmatrix((double *)R_alloc(nstate*nstate, sizeof(double)),
                               nstate, nstate);
    chaz = (double**) dmatrix((double *)R_alloc(nstate*nstate, sizeof(double)),
                               nstate, nstate);
    if (sefit >0)  
        umat = (double**) dmatrix((double *)R_alloc(nperson*nstate, sizeof(double)),
                               nstate, nperson);

    /* R_alloc does not zero allocated memory */
    for (i=0; i<nstate; i++) {
        ws[i] =0;
        ns[i] =0;
        nev[i] =0;
        for (j=0; j<nstate; j++) {
                hmat[i][j] =0;
                chaz[i][j] =0;
        }
    }
    for (i=0; i<nperson; i++) {
        atrisk[i] =0;
        wtp[i] = 0.0;
        dstate[i] = cstate[i];  /* cstate starts as the initial state */
    }
    if (sefit ==1 || influence >0) {
        dptr = i0;
        for (j=0; j<nstate; j++) {
            for (i=0; i<nperson; i++) umat[i][j] = *dptr++;
        }
     }
    if (influence & 0x1) {
        /* The initial influence is saved out first */
        dptr = i0;
        for (j=0; j<nstate; j++) {
            for (i=0; i<nperson; i++) 
                *imat1++ = *dptr++;   /* save in the output */
        }
    } 
    if (influence & 0x4) {
        /* Zero the first time point, since we accumulate results */
        for (j=0; j<(nstate*nperson); j++) {
            imat3[j] = 0;
        }
    } 
    itime =0; /*current time index, for output arrays */
    eptr  = 0; /*index to sort1, the entry times */
    itimej= 0; /* current index for the inftimes vector */
    ctime = etime[sort2[0]];
    for (i=0; i<n; ) {
        ck = sort2[i];
        ctime = etime[ck];  /* current time value of interest */
        if (influence) {
            for (; itimej < ninftime && inftime[itimej] < ctime; itimej++) {
                if (influence & 0x1) {
                    for (j=0; j<nstate; j++) {
                        for (k=0; k<nperson; k++) *imat1++ = umat[k][j];
                    }
                }
                if (influence & 0x4) {
                    if (inftime[itimej] > lagtime) {
                        /* this chunk only runs inftimes that are between event times */
                        for (j=0; j<nstate; j++) {
                            for (k=0; k<nperson; k++) imat3[k + j*nperson] += umat[k][j] *
                                                          (inftime[itimej] - lagtime);
                        }
                        lagtime = inftime[itimej];
                    }
                    if ((itimej+1) < ninftime) {
                        /* step ahead */
                        k = nstate* nperson;
                        for (j=0; j< k; j++)
                            imat3[k+j] = imat3[j];
                        imat3 += k;
                    }
                }
            }
            if ((influence & 0x4) && (ctime > lagtime)) {
                /* the inf1 value is about to change, catch up */
                for (j=0; j<nstate; j++) {
                    for (k=0; k<nperson; k++) imat3[j*nperson +k] += umat[k][j] * 
                                                      (ctime - lagtime);
                }
                lagtime = ctime;
            }
        }

        /* Add subjects whose entry time is < ctime into the counts */
        for (; eptr<n; eptr++) {
            k = sort1[eptr];
            if (entry[k] < ctime) {
                kk = cstate[id[k]];  /*current state of the addition */
                ns[kk]++;
                ws[kk] += wt[k];
                wtp[id[k]] = wt[k];
                atrisk[id[k]] =1;   /* mark them as being at risk */
            }
            else break;
        }
            
        for (j=0; j<nstate; j++) {
            for (k=0; k<nstate; k++) {
                hmat[j][k] =0;
            }
         }

        /* Count up the number of events and censored at this time point */
        for (k=0; k<nstate; k++) nev[k] =0;
        ncensor[itime] =0;
        wevent =0;
        for (j=i; j<n; j++) {
            k = sort2[j];
            if (etime[k] == ctime) {
                if (status[k] >0) {
                    newstate = status[k] -1;  /* 0 based subscripts */
                    oldstate = cstate[id[k]];
                    if (oldstate != newstate) {
                        /* A "move" to the same state does not count */
                        dstate[id[k]] = newstate;
                        nev[newstate]++;
                        wevent += wt[k];
                        hmat[oldstate][newstate] += wt[k];
                    }
                }
                else ncensor[itime]++;
            }
            else break;
         }
                
        if (wevent > 0) {  /* there was at least one move with weight > 0 */
            /* finish computing H */
            for (j=0; j<nstate; j++) {
                if (ns[j] >0) {
                    temp =0;
                    for (k=0; k<nstate; k++) {
                        temp += hmat[j][k];
                        hmat[j][k] /= ws[j];  /* events/n */
                    }
                    hmat[j][j] =1 -temp/ws[j]; /*rows sum to one */
                }
                else hmat[j][j] =1.0; 
         
            }
            if ((sefit >0) || (influence >0)) {
                if (sefit >0) {
                    /* Update U, part 1  U = U %*% H -- matrix multiplication */
                    for (j=0; j<nperson; j++) { /* row of U */
                        for (k=0; k<nstate; k++) { /* column of U */
                            temp2[k]=0;
                            for (kk=0; kk<nstate; kk++) 
                                temp2[k] += umat[j][kk] * hmat[kk][k];
                        }  
                        for (k=0; k<nstate; k++) umat[j][k] = temp2[k];
                    }

                    /* step 2, add in dH term */
                    for (j=0; j<nperson; j++) {
                        if (atrisk[j]==1) {
                            oldstate = cstate[j];
                            for (k=0; k<nstate; k++)
                                umat[j][k] -= hmat[oldstate][k]* p[oldstate]/ ws[oldstate];
                            umat[j][dstate[j]] += p[oldstate]/ws[oldstate];
                        }
                    }
                }
            }
            /* Finally, update chaz and p.  */
            for (j=0; j<nstate; j++) {
                for (k=0; k<nstate; k++) chaz[j][k] += hmat[j][k];
                chaz[j][j] -=1;  /* Update using H2 */

                temp2[j] =0;
                for (k=0; k<nstate; k++)
                    temp2[j] += p[k] * hmat[k][j];
             }
            for (j=0; j<nstate; j++) p[j] = temp2[j];
        }
        /* store into the matrices that will be passed back */
        for (j=0; j<nstate; j++) {
            pmat[j][itime] = p[j];
            nrisk[j][itime] = ns[j];
            nevent[j][itime] = nev[j];
            for (k=0; k<nstate; k++) *cumhaz++ = chaz[k][j];
            if (sefit >0) {
                temp =0;
                for (k=0; k<nperson; k++) 
                    temp += wtp[k]* wtp[k]*umat[k][j]*umat[k][j];
                vmat[j][itime] = sqrt(temp);
            }
        }
      
        /* Take the current events and censors out of the risk set */
        for (; i<n; i++) {
            j= sort2[i];
            if (etime[j] == ctime) {
                oldstate = cstate[id[j]]; /*current state */
                ns[oldstate]--;
                ws[oldstate] -= wt[j];
                if (status[j] >0) cstate[id[j]] = status[j]-1; /*new state */
                atrisk[id[j]] =0;
            }
            else break;
        }
        itime++;  
    }

    if (influence) {
        for (; itimej < ninftime; itimej++) {
            if (influence & 0x1) {
                for (j=0; j<nstate; j++) {
                    for (k=0; k<nperson; k++) *imat1++ = umat[k][j];
                }
            }
            if (influence & 0x4) {
                if (inftime[itimej] > lagtime) {
                    for (j=0; j<nstate; j++) {
                        for (k=0; k<nperson; k++) imat3[k + j*nperson] += umat[k][j] *
                                                      (inftime[itimej] - lagtime);
                    }
                    lagtime = inftime[itimej];
                }
                if ((itimej+1) < ninftime) {
                    /* step ahead */
                    k = nstate* nperson;
                    for (j=0; j< k; j++)
                        imat3[k+j] = imat3[j];
                    imat3 += k;
                }
            }
        }
    }
    /* return a list */
    UNPROTECT(3);
    return(rlist);
}
