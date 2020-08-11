/*
** This is a modified copy of survfitci, whose sole job is to compute
**  the influence matrix of the AJ estimate, at a particular set of points
** Extra arguments are the number and location of reporting times, and
**  whether mean time in state is required.
**
** Arguments:  For an R objects 'zed' I use 'zed2' to refer to the object
**     and 'zed' for the contents of the object.
**  ftime = two column matrix containing entry and exit
**  status= vector of status codes, at end of each interval
**  sort1 = order vector for the entry times
**  sort2 = order vector for the exit times
**  cstate= current state for each observation
**  wt    = case weight for each observation
**  curve = which curve for the obs.  Data is sorted by time within curve,
**   and accumulation starts afresh at a new one.
**  p0    = initial distribution of states
**  i0    = initial influence matrix, will be number of subjects by number
**           of states, often will be 0.
**  otime = vector of output times
**  starttime = starting time for the curve, needed for AUC
**  doauc = 1 if AUC is needed, 0 otherwise
*/
#include "survS.h"
#include "survproto.h"
#include <math.h>

SEXP survfitresid(SEXP ftime2,  SEXP status2, SEXP sort12,     SEXP sort22, 
	          SEXP cstate2, SEXP wt2,     SEXP curve2,     SEXP p02,     
	          SEXP i02,     SEXP otime2,  SEXP starttime2, SEXP doauc2) {   
    int i, j, k, kk;   /* generic loop indices */
    int p1, p2;        /* observations pointed to by sort1 and sort2 */
    int nout;          /* number of output times */
    int itime, eptr, psave; /*specific indices */
    int icurve;        /* the current curve being processed */
    double ctime;      /*current time of interest, in the main loop */
    int oldstate, newstate; /*when changing state */

    double temp, temp2, *tempvec;  /* scratch doubles and vector */
    double *dptr;      /* reused in multiple contexts */
    double *p;         /* current prevalence vector */
    double **hmat;      /* hazard matrix at this time point */
    double **cmat;     /* H = I + C, this is C */

    int *atrisk;       /* 1 if the observation is currently at risk */
    int   nevent;      /* number of events at this time */
    double *ws;        /* weighted count at risk, by state */
    double wevent;     /* weighted number of events at current time */
    int nstate;        /* number of states */
    int nobs;          /*number of obs */
    double starttime;
    int doauc;

    /* pointers to the R variables */
    int *sort1, *sort2;  /*sort index for entry time, event time */
    double *entry,* etime;  /*entry time, event time */
    int *status;        /*0=censored, 1,2,... new states */
    int *cstate;        /* current state for each subject */
    double *wt;         /* weight for each observation */
    double *i0;         /* initial influence */
    int *curve;         /* for each obs, which curve is it in */
    int *otime;         /* requested output times */
    SEXP setemp;        /* a temp for variables on the list */
        
    /* returned objects */
    SEXP rlist;         /* the returned list and variable names of same */  
    const char *rnames[]= {"influence.pstate", "influence.auc", ""}; 
    double **infa, **infp; /* pointers to influence arrays */

    nobs   = LENGTH(sort12);    /* number of observations in the data */
    cstate  = INTEGER(cstate2);
    entry= REAL(ftime2);
    etime= entry + nobs;        /* don't overlap name with exit() */
    status= INTEGER(status2);
    sort1= INTEGER(sort12);
    sort2= INTEGER(sort22);
    wt = REAL(wt2);
    curve = INTEGER(curve2);
    nstate = LENGTH(p02);  /* number of states */
    i0 = REAL(i02);
    otime = INTEGER(otime2);
    nout  = LENGTH(otime2);  /* number of output times */
    starttime = REAL(starttime2)[0]; /* single number vector */
    doauc = asInteger(doauc2);

    /* allocate space for the output objects
    ** Ones that are put into a list do not need to be protected
    */
    PROTECT(rlist=mkNamed(VECSXP, rnames));
    setemp = SET_VECTOR_ELT(rlist, 0, allocMatrix(REALSXP, nobs, nout*nstate));
    infp   = dmatrix(REAL(setemp), nobs, nstate);
    if (doauc==1) {
	setemp = SET_VECTOR_ELT(rlist, 1, allocMatrix(REALSXP, nobs, nout*nstate));
	infa   = dmatrix(REAL(setemp), nobs, nstate); /* influence for AUC */
    }
    
    /* allocate space for scratch vectors */
    ws = (double *) R_alloc(3*nstate, sizeof(double)); /*weighted number in state */
    p = ws + nstate;
    tempvec = p + nstate;

    atrisk = (int *) R_alloc(nobs, sizeof(int));
    hmat = (double**) dmatrix((double *)R_alloc(nstate*nstate, sizeof(double)),
                               nstate, nstate);
    cmat = (double**) dmatrix((double *)R_alloc(nstate*nstate, sizeof(double)),
                               nstate, nstate);

    /* R_alloc does not zero allocated memory */
    for (i=0; i<nstate; i++) {
        ws[i] =0;
        for (j=0; j<nstate; j++) {
                hmat[i][j] =0;
                cmat[i][j] =0;
        }
    }
    for (i=0; i<nobs; i++) atrisk[i] =0;
    dptr = i0;
    for (j=0; j<nstate; j++) {
	for (i=0; i<nobs; i++) infp[j][i] = *dptr++;
    } 
    
    eptr  = 0; /*index to sort1, the entry times */
    i = 0;
    while (i<nobs) {
	icurve = curve[sort2[i]];
	itime  =0;
	for (; i<nobs; i++) {
	    p2 = sort2[i];
	    if (curve[p2] != icurve) break;
	    ctime = etime[p2];  /* current time value of interest */
	
		/* this is cast as a loop in case there are multipe output times
		**  between two data times
		*/
		for(; otime[itime] <= ctime && (itime < nout); itime++) {
		    /* finished with one of the output times, update */
		    if (doauc ==1) {
			/* AUC influence += (Phat influence)*(ctime- starttime) */
			for (j=0; j<nobs; j++) {
			    for (k=0; k<nstate; k++)
				infa[k][j] += infp[k][j]* (otime[itime]-starttime);
			}
			starttime = otime[itime];
		    }

		    /* Copy it forward to the next output time 
		    ** The matrix indices infa and infp point to the n x nstate
		    **  influences for the current output time.  We copy forward
		    **  to the next array slice and update the pointers
		    */
		    if ((itime +1) < nout) {
			for (k=0; k<nstate; k++){
			    tempvec = infp[k];
			    infp[k] += nstate*nobs;
			    for (j=0; j<nobs; j++)
				infp[k][j] = tempvec[j];
			}
			if (doauc ==1) {
			    for (k=0; k<nstate; k++){
				tempvec = infa[k];
				infa[k] += nstate*nobs;
				for (j=0; j<nobs; j++)
				    infa[k][j] = tempvec[j];
			    }	
			}
		    }
		}
		if (itime== nout) break;  /* no need to go past last output time*/

		/* Add subjects whose entry time is < ctime into the counts */
		for (; eptr<nobs; eptr++) {
		    p1 = sort1[eptr];
		    if (entry[p1] < ctime) {
			k = cstate[p1];  /*entry state of the addition */
			ws[k] += wt[k]; /* weighted number at risk */
			atrisk[p1] =1;   /* mark them as being at risk */
		    }
		    else break;
		}	
            
		for (j=0; j<nstate; j++) {
		    for (k=0; k<nstate; k++) {
			cmat[j][k] =0;
		    }
		}
       
		/* Count up the number of transitions at this time point */
		nevent=0; 
		wevent =0;
		for (j=i; j<nobs; j++) {
		    p2 = sort2[j];
		    if (etime[p2] > ctime) break;
		    if (status[p2] >0 && (cstate[p2] != status[p2])) {
			/* A "move" to the same state does not count */
			newstate = status[p2] -1;  /* 0 based subscripts */
			oldstate = cstate[p2];
			psave = p2;
			nevent++;
			wevent += wt[p2];
			cmat[oldstate][newstate] += wt[p2]/ ws[oldstate];
			cmat[oldstate][oldstate] -= wt[p2]/ ws[oldstate];
		    }
		}
             
		if (nevent >0 && doauc==1) {
		    /* AUC influence uses the prior infa matrix, 
		       before we update infa, but
		       only needs to be updated at an event */
		    temp = ctime -starttime;
		    for (j=0; j< nobs * nstate; j++) infa[itime][j] += infp[itime][j]*temp;
		    starttime = ctime;
		}

		if (nevent ==1) {
		    /*
		    ** The H matrix is actually I + C, we have so far computed C
		    ** In this single event case, which is moderately common,
		    **  we can use a faster update due to the simplicity of C
		    */
		    temp = p[oldstate]*cmat[oldstate][oldstate];
 
		    for (j=0; j<nobs; j++) {
			/* update U, part 1, new U = UH = U + UC  (U = infa)*/
			infp[oldstate][i] += temp* infp[oldstate][i];
			infp[newstate][i] -= temp* infp[newstate][i];
		    }
		
		    /* add partial C/wt[i], which affects all those in oldstate*/
		    temp2 = p[oldstate]/ws[oldstate];
		    infp[oldstate][psave] += temp2;
		    infp[newstate][psave] -= temp2;
		    for (j=i; j<nobs; j++) {
			p2 = sort2[i];
			if (atrisk[p2] && cstate[p2] == oldstate) {
			    infp[oldstate][p2] -= temp/ws[oldstate];
			    infp[newstate][p2] += temp/ws[oldstate];
			}
		    }
		} else if (nevent >1) {
		    /* Update U, part 1  U = U + U %*% C -- matrix multiplication */
		    for (j=0; j<nobs; j++) { /* row of U */
			for (k=0; k<nstate; k++) { /* column of U */
			    tempvec[k]=0;
			    for (kk=0; kk<nstate; kk++) 
				tempvec[k] += infp[j][kk] * cmat[kk][k];
			}  
			for (k=0; k<nstate; k++) infp[j][k] += tempvec[k];
		    }

		    /* step 2, add in dH term */
		    for (j=i; j<nobs; j++) {
			p2 = sort2[j];
			if (atrisk[p2]==1) {
			    oldstate = cstate[p2];
			    for (k=0; k<nstate; k++)
				infa[k][p2] -= cmat[oldstate][k]* p[oldstate]/ ws[oldstate];
			    if (oldstate != status[p2]) { /* transition*/
				infa[oldstate][p2] += p[oldstate]/ws[oldstate];
				infa[status[p2]][p2]-= p[oldstate]/ws[oldstate];
			    }
			}
		    }

		    /* Update p  */
		    for (j=0; j<nstate; j++) {
			tempvec[j] =0;
			for (k=0; k<nstate; k++)
			    tempvec[j] += p[k] * cmat[k][j];
		    }
		    for (j=0; j<nstate; j++) p[j] += tempvec[j];
		}
 
		/* Take the current events and censors out of the risk set */
		for (; i<nobs; i++) {
		    p2 = sort2[i];
		    if (etime[p2] == ctime) {
			oldstate = cstate[p2]; /*current state */
			ws[oldstate] -= wt[p2];
			atrisk[p2] =0;
		    }
		    else break;
		}
	    } 

	for (; itime < nout; itime++) {
	    /* there are reporting times after the last event; finish up */
	    if (doauc==1) {
		for (k=0; k<nstate; k++) {
		    for (j=0; j<nobs; j++) 
			infa[k][j] += infp[k][j] * (otime[itime]-starttime);
		}
		starttime = otime[itime];
	    }

	    if ((itime+1)<nout) {
		/* there is yet one more reporting time after this one */
		for (k=0; k<nstate; k++){
		    tempvec = infp[k];
		    infp[k] += nstate*nobs;
		    for (j=0; j<nobs; j++)
			infp[k][j] = tempvec[j];
		}	
		if (doauc ==1) {
		    for (k=0; k<nstate; k++){
			tempvec = infa[k];
			infa[k] += nstate*nobs;
			for (j=0; j<nobs; j++)
			    infa[k][j] = tempvec[j];
		    }	
		}
	    }
	}
    }	
    UNPROTECT(1);
    return(rlist);
}
