/*
** Code for the Aalen-Johansen estimate, and influence matrices for the same
**
**  I tend to use "zed2" for the R variable, and then "zed" for the C pointer
**  to its contents
** y     : matrix of survival, columns of time1, time2, state
** sort1 : sort vector for time1
** sort2 : sort vector for time2
** utime : vector of output times; estimates are created for these times
** cstate: state during each interval, before any transition at the end
** wt    : case weight
** grp   : group, for the grouped variance (very often the subject id)
**        this is always 0, 1, 2, ...
** ngrp  : number of unique groups
** position : when a subject has multiple obs in sequence (t1, t2), (t2, t3),
**         ..., this is 1 for the first, 2 for the last, 3 for one that is
**          both first and last. This drives the entry and censoring counts
** p0    : initial vector p(0)
** i0    : initial influence matrix
** sefit : 0= no se, 1=se, 2= influence matrix returned
** entry : if 1, then utime contains entry times, and we return entry counts
** hindx : matrix that indexes the hazards: the [k,j] element is the output
**         column for the j:k transtion  (given j,k; find the col of chaz)
** trmat : two column matrix with one row per transition, with elements of 
**          from:to (given the col of chaz, find j and k)
** t0    : the starting point for the curve, needed for std of the AUC
**
** When there are multiple curves, e.g., survfit(Surv(t1, t2, state) ~ z) this
**  routine is called once for each value of z
** y, wt, cstate, grp, pattern will have values for all subjects
** sort1, sort2 will have nperson rows: all curve 2 for instance, and
**  point to the rows of interest in y
** ngrp will be the number of unique ids in the subset
*/

#include "survS.h"
#include "survproto.h"
#include <math.h>

SEXP survfitaj(SEXP y2,       SEXP sort12,  SEXP sort22, SEXP utime2,
               SEXP cstate2,  SEXP wt2,     SEXP grp2,   SEXP ngrp2,
	       SEXP p02,      SEXP i02,     SEXP sefit2, SEXP entry2,
	       SEXP position2,SEXP hindx2,  SEXP trmat2, SEXP t02) {   
    int i, j, k, jk, g, p;   /* generic loop indices */
    int i1, i2, person1, person2;  /* index the sort1/sort2 vectors */
    double ctime;      /*current time of interest, in the main loop */
    int tdeath;        /* deaths at a given time point */
    int nstate;        /* number of states */
    int n, nused;     /* number of rows in y, rows in sort1/sort2 */
    int ngrp;          /* number of groups (clusters) */
    double delta;      /* temp, width of time interval */
    double temp, haz, htemp;

    /* working variables, memory will be obtained using Ralloc */
    /* All the "=0" ones are allocated/used only if the robust variance is
       needed, and gcc warns "may be unitialized" without the =0 
       The warnings are all false positives, but CRAN will be unhappy
    */
    double **U=0, **UA=0;  /* dfbeta matrices for pstate and AUC */
    double **C=0;          /* dfbeta matrix for cumhaz, one row per group */
    double **H=0;          /* H matrix (nstate x nstate) */
    double **Ucopy=0;      /* temp copy of U, for matrix multiplication */
    double *hcopy;         /* bookkeeping for matrix mult */
    double *se1=0, *se2=0, *se3=0;  /* sums of squares, for variances */
    double **wg=0;        /* weighted number at risk, by cluster and state */
    double *ntemp;        /* running total for number at risk, per state */
    double *phat;
    double *chaz;

    /* pointers to the R input variables */
    int *sort1, *sort2;  /*sort indexs for time1 and time2 */
    double *time1, 
	   *time2,
	   *state;      /* (time1, time2, state) survival time */
    int ntime;          /* number of unique event time values */
    int *cstate;        /* current state for each subject */
    double *wt;         /* weight for each observation */
    double *p0;         /* initial state probabilities */
    double *i0;         /* initial influence */
    int *grp;           /* for each obs, which group for variance */
    int sefit;          /* 0 no std err, 1=stderr, 2= std and influence.pstate*/
    int entry, nhaz;    /* flag: do entry counts?   Number of hazards */
    int *hindx, *trmat;
    int *position;       /* identify un-informative breaks */
    double *utime;      /* reporting times for the result */
    double t0;          /* starting point for AUC */
        
    /* returned objects for R */
    SEXP rlist;         /* the returned list and variable names of same */  
    const char *rnames[]= {"n.risk", "n.event", "n.censor", "pstate", "cumhaz", 
	                   "std.err", "std.chaz", "std.auc", "influence", 
			   "n.enter", "n.transition", ""};
    SEXP setemp;
    double **pstate, **cumhaz, *usave=0; /*=0 to silence -Wall warning */
    double **nrisk, **nevent, **ncensor, **nenter=0, **ntrans;
    double **stdp,  **stdc, **stda; 

    ntime= LENGTH(utime2);
    n    = nrows(y2); /* number of unique subjects */
    nused = LENGTH(sort12);    /* number of observations used for this curve */
    time1= REAL(y2);
    time2= time1 + n;
    state = time2 + n;  /* the last use of 'n' in this program */
    sort1= INTEGER(sort12);
    sort2= INTEGER(sort22);
    utime= REAL(utime2);
    cstate = INTEGER(cstate2);
    wt = REAL(wt2);
    ngrp = asInteger(ngrp2);
    grp  = INTEGER(grp2);
    p0 = REAL(p02);
    nstate = LENGTH(p02);  /* number of states */
    i0 = REAL(i02);
    sefit = asInteger(sefit2);
    entry = asInteger(entry2);
    nstate = nrows(hindx2); /* number of states */
    nhaz   = nrows(trmat2);  /* number of unique transtion types*/
    hindx =  INTEGER(hindx2); /* the [from, to] element contains trans index*/
    trmat  =  INTEGER(trmat2);  /* row k is the from:to for transition k */
    position= INTEGER(position2);
    t0 = asReal(t02);

    /*
    ** If we know from:to and need the transition number, 
    **   hindx[from + nstate *to] contains it
    ** If we have the transition number and need from:to, trmat[transition]
    **   has 'from', and trmat[transtion + ntran] has 'to'
    ** A model with 5 states has 25 potential transitions, but only a small
    **   handful of those may actually occur. nstate and ntrans.
    */
 
    /* allocate space for the output objects
    ** Objects that are put into a list do not need to be protected, only
    **  the list itself need be.
    ** I return both weighted and unweighted counts.  The weighted first in a
    **  vector, then unweighted.
    */
    PROTECT(rlist=mkNamed(VECSXP, rnames));
    /* n.risk = counts for each state */
    setemp = SET_VECTOR_ELT(rlist, 0, allocMatrix(REALSXP, ntime, 
						  2*nstate));
    nrisk = dmatrix(REAL(setemp), ntime, 2*nstate);
    /* n.transition has a column per transition, n.event one per state*/
    setemp = SET_VECTOR_ELT(rlist, 1, allocMatrix(REALSXP, ntime, nstate));
    nevent = dmatrix(REAL(setemp), ntime, nstate);
    setemp = SET_VECTOR_ELT(rlist, 10, allocMatrix(REALSXP, ntime, 2*nhaz));
    ntrans = dmatrix(REAL(setemp), ntime, 2*nhaz);
    /* censor is per state, ditto for entry */
    setemp  = SET_VECTOR_ELT(rlist, 2, allocMatrix(REALSXP, ntime, 2*nstate));
    ncensor = dmatrix(REAL(setemp), ntime, 2*nstate);
    if (entry==1) {
	setemp = SET_VECTOR_ELT(rlist,9, allocMatrix(REALSXP, ntime, 2*nstate)); 
	nenter = dmatrix(REAL(setemp), ntime, 2*nstate);
    }
    setemp = SET_VECTOR_ELT(rlist, 3, allocMatrix(REALSXP, ntime, nstate));
    pstate = dmatrix(REAL(setemp), ntime, nstate);
    setemp = SET_VECTOR_ELT(rlist, 4, allocMatrix(REALSXP, ntime, nhaz));
    cumhaz = dmatrix(REAL(setemp), ntime, nhaz);
    /* memory must be zeroed */
    for (j=0; j< nstate; j++) {
	for (i=0; i<ntime; i++) {
	    nrisk[j][i] =0;
	    nrisk[j+nstate][i] =0;
	    ncensor[j][i] =0; 
	    ncensor[j+nstate][i] =0;
	    nevent[j][i] =0;
	}
	if (entry==1) {
	    for (i=0; i<ntime; i++)  {
		nenter[j][i] =0;
		nenter[j+nstate][i] =0;
	    }
	}
    }	
    for (jk=0; jk < nhaz; jk++) {
	for (i=0; i<ntime; i++) {
	    ntrans[jk][i] =0;
	    ntrans[jk+ nhaz][i] =0; 
	    cumhaz[jk][i] =0;
	}
    }

    if (sefit >0) {
        setemp = SET_VECTOR_ELT(rlist, 5,  allocMatrix(REALSXP, ntime, 
						       nstate));
        stdp = dmatrix(REAL(setemp), ntime, nstate); 
        setemp = SET_VECTOR_ELT(rlist, 6,  allocMatrix(REALSXP, ntime,
						       nhaz));
	stdc = dmatrix(REAL(setemp), ntime, nhaz); /*cumhaz */
        setemp = SET_VECTOR_ELT(rlist, 7,  allocMatrix(REALSXP, ntime, 
						       nstate));
	stda = dmatrix(REAL(setemp), ntime, nstate); /* AUC */
    }
    /* 
    ** If needed, allocate space for a returned influence matrix.
    ** They can be huge.
    */
    if (sefit>1) {
	/* Allocate the R object as a matrix, which can be larger than a vector:
	 nrow, ncol, and vector lengths are all integers */
        setemp =SET_VECTOR_ELT(rlist,8,allocMatrix(REALSXP, ngrp*nstate,ntime));
	usave = REAL(setemp);
	/* usave is a pointer to a double, which is bigger than integer and
	   won't overflow as we increment it.
	   All elements will be written, so no need to zero it.
	*/
    }

    /* allocate space for scratch vectors */
    H =  dmatrix((double *)  R_alloc(nstate*nstate, sizeof(double)), nstate,
		 nstate);
    ntemp = (double *) R_alloc(4*nstate, sizeof(double));
    phat  = ntemp + 2*nstate;
    hcopy = phat + nstate;
    chaz  = (double *) R_alloc(nhaz, sizeof(double));

    for (i=0; i<4*nstate; i++) ntemp[i] =0;
    for (i=0; i< nhaz; i++) chaz[i]=0;

    if (sefit > 0) {
	U  = dmatrix((double *) R_alloc(ngrp*nstate, sizeof(double)), ngrp, 
		     nstate);
	Ucopy  = dmatrix((double *) R_alloc(ngrp*nstate, sizeof(double)), ngrp, 
			 nstate);
	UA = dmatrix((double *) R_alloc(ngrp*nstate, sizeof(double)), ngrp, 
		     nstate);
	C  = dmatrix((double *) R_alloc(ngrp*nhaz,   sizeof(double)), ngrp,
		     nhaz);
	wg = dmatrix((double *) R_alloc(ngrp*nstate, sizeof(double)), ngrp, 
		     nstate);
	se1 = (double *) R_alloc(2*nstate + nhaz, sizeof(double));
	se2 = se1 + nstate;
	se3 = se2 + nhaz;

	/* R_alloc does not zero allocated memory */
	for (i=0; i<nstate; i++) {
	    se1[i] =0;
	    se2[i] =0;
	    se3[i] =0;
	    for (g=0; g<ngrp; g++) {
		U[i][g] = 0;
		UA[i][g] =0;
		wg[i][g] =0;
	    }
	}
	for (jk=0; jk<nhaz; jk++) { 
	    for (g=0; g<ngrp; g++) C[jk][g] =0;
	}
    }

    /*
    ** First pass through the data, to create all the counts.
    **  This is done in reverse time order: for most data, risk sets
    **  decrease over time, so this direction is a bit more accurate (more 
    **  additions than subtractions).   
    ** ntemp is the number at risk, weighted and unweighted (so of length
    **  2*nstate) and is accumulated. It is driven by the cstate 
    **  variable, i.e. current state in each interval (observation).
    ** Other counts are at this time point and do not accumulate.
    ** 	
    ** person1 tracks sort1, person2 tracks sort2
    ** Assume a subject with obs of (2, 5a],(5,6+], (6,8b], [8, 10a], (10,11+]
    **   cstate of s0, a, a, b, a
    ** At time 11 add one to state a, at 10 remove from a and add to b,
    **  at 8 remove from b and add to a, at 6 subtract weight from a and then
    **  add to a (the weight might change), at 5 remove from a and
    **  add to s0, at 2 remove from s0. 
    ** See methods, survfit:number at risk.
    */
    for (j=0; j< 2*nstate; j++) ntemp[j]=0;
    person1 = nused-1; person2= nused -1;
    for (i=ntime-1; i>=0; i--) {
	ctime = utime[i];  /* current time point of interest */
	/*
	** Utime moves from right (later times) to earlier, with ctime the
	**  most recent checkpoint:
	** First remove weights for intervals that no longer overlap ctime,
	**  i.e., time1 >= ctime.  If entry=1 and this is the first interval
	**  for a subject (position =1 or 3) add to the entry count.  Note that
	**  if entry >0, that all entry times will be in the utime vector.
        */
	for (; person1 >=0 && time1[sort1[person1]] >= ctime; person1--) {
	    i1 = sort1[person1];
	    j = cstate[i1];
	    ntemp[j] -= wt[i1];
	    ntemp[j+nstate]--;
	    if (position[i1] & 01 & entry) { /* start of a sequence */
		nenter[j][i] += wt[i1];
		nenter[j+nstate][i]++;
	    }
	}
	/* 
	** Now deal with intervals that overlap ctime, but did not at prior iter,
	**   time2 >= ctime, and add them to the risk set. 
	** If it is an event add it to nevent, if it is a non-event and
	**   position = 2 or 3, then add it to the censoring count.
	*/  
	for (; person2>=0 && time2[sort2[person2]] >= ctime; person2--) {
	    i2 = sort2[person2];
	    j = cstate[i2];  /* current state*/
	    ntemp[j] += wt[i2];
	    ntemp[j+nstate]++;
	    if (state[i2] >0) { /* an event*/
		k = state[i2] -1;  /* the event type */
		jk = hindx[ j + k*nstate]; /* look up transition number */
		ntrans[jk][i] += wt[i2];
		ntrans[jk+ nhaz][i]++;
		nevent[k][i] += wt[i2];
	    }
	    else if (position[i2] >1 ){ 
		/* censored, and last interval for this subject*/
		ncensor[j][i] += wt[i2];
		ncensor[j+nstate][i]++;
	    } 
	}
	for (j=0; j< 2*nstate; j++) nrisk[j][i] = ntemp[j]; /*save result*/
    }

    /* Initial estimates.  UA and C start at 0, initial for U was passed in*/
    for (j=0; j<nstate; j++) phat[j] = p0[j];
    if (sefit >0) {
	for (j=0; j<nstate; j++) {
	    temp =0;
	    for (g=0; g< ngrp; g++) {
	        U[j][g] = i0[g + j*ngrp];
		wg[j][g] =0;	/* weighted number at risk, by group */		
		temp += U[j][g]* U[j][g];
  	    }
	    se1[j] = sqrt(temp);
	}
    }

    /*
    ** Now walk forward in time and compute the AJ
    ** Formula documentation is in the design and methods document,
    **  see Andersen-Gill: influence.
    **  Each transition jk from state j to state k causes an increment
    **  jk = column of the chaz and ntrans matrices
    */
    person1 =0;
    person2 =0;
    for (i=0; i<ntime; i++){ 
	if (sefit >0) {
	    /* update the AUC */
	    if (i>0) delta = utime[i] - utime[i-1];
	    else delta = utime[i] - t0;
	    for (j=0; j< nstate; j++) {
		temp =0;
		for (g=0; g<ngrp; g++) {
		    UA[j][g] += delta*U[j][g];
		    temp += UA[j][g] * UA[j][g];
		}
		se3[j]= sqrt(temp);
	    }
	    /* 
	    ** update wg = sum of those at risk, at this time, in each group
	    ** At the end
	    **   sort1[person1] will point to the rightmost start time < utime[i]
	    **   sort2[person2] to the leftmost ending time that is < utime[i]
	    */
	    for (; person1 <nused; person1++) {
		i1= sort1[person1];
		if (time1[i1] < utime[i]) wg[cstate[i1]][grp[i1]] += wt[i1];
		else break;
	    }
	    for(;  person2 <nused; person2++) { 
		i2 = sort2[person2];
		if (time2[i2] < utime[i]) wg[cstate[i2]][grp[i2]] -= wt[i2];
		else break;
	    }

	    /* Update the dN part of C, and compute H (acually H-I) */
	    tdeath =0;
	    for (j=0; j<nstate; j++) {
		hcopy[j] =0;
		for (k=0; k<nstate; k++) H[j][k] =0;
	    }
	    for (p=person2; p < nused; p++) {
		i2 = sort2[p];
		if (time2[i2] > utime[i]) break;  /* later */
		if (time2[i2] == utime[i] && state[i2] > 0) {
		    tdeath++;
		    j= cstate[i2];
		    k= state[i2] -1;  /* first level is 'censored' */
		    jk = hindx[j + k*nstate];		    
		    g= grp[i2];	
		    C[jk][g] += wt[i2]/nrisk[j][i];
		    if (j != k) {
			H[j][j] -= wt[i2]/nrisk[j][i];
			H[j][k] += wt[i2]/nrisk[j][i];
		    }
		}
	    }
	    if (tdeath==0) goto saveit; /* no events, C and U won't change*/
	    
	    /* 
	    ** set U = U+ UH, matrix multiplication, take advantage that H is
	    **  normally very sparse, often only 2 non-zero elements. 
	    ** Off diagonal elements are 0 or positive, H[j,j]==0 implies the
	    **  entire row is zero.
	    ** Algorithm:
	    **   If H[j,k] !=0, then increment U[k,] += U[j,]* H[j,k]
	    **   But before modifying U[k,] look ahead to see if we we will need
	    **    an unmodified copy later as the "U[j,]" part of the above
	    **    equation, i.e., if H[k,k] !=0.   Only copy when we must.
	    ** 
	    */
	    for (j=0; j<nstate; j++) {
		if (H[j][j] !=0) { /* a non-zero row */
		    for (k=0; k<nstate; k++) {
			if (k!=j && H[j][k]!=0) {
			    if (H[k][k] !=0 && k>j && hcopy[k] ==0) {
				/* save a copy of this column */
				hcopy[k] =1;
				for (g=0; g<ngrp; g++)
				    Ucopy[k][g] = U[k][g];
				}
			    if (hcopy[j] == 0) {
				for (g=0; g<ngrp; g++) 
				    U[k][g] += U[j][g]* H[j][k];
			    } else {
				for (g=0; g<ngrp; g++) 
				    U[k][g] += Ucopy[j][g]* H[j][k];
			    }
			}
		    }
		    if (hcopy[j] ==0)
			for (g=0; g<ngrp; g++) U[j][g] += U[j][g]*H[j][j];
		    else
			for (g=0; g<ngrp; g++) U[j][g] += Ucopy[j][g]*H[j][j];
		}
	    }

	    /* Update the dN part of U */
	    for (p=person2; p < nused; p++) {
		i2 = sort2[p];
		if (time2[i2] > utime[i]) break;  /* later */
		if (time2[i2] == utime[i] && state[i2] > 0) {
		    j= cstate[i2];
		    k= state[i2] -1;  /* first level is 'censored' */
		    g= grp[i2];	
		    if (j != k) {
			U[j][g] -= wt[i2]*phat[j]/nrisk[j][i];
			U[k][g] += wt[i2]*phat[j]/nrisk[j][i];
		    }
		}
	    }

	    /* Update the phat part of the influence */
	    for (jk=0; jk<nhaz; jk++) {
		if (ntrans[jk][i] > 0) {
		    /* make updates to the working estimates */
		    j = trmat[jk]; k = trmat[jk + nhaz];  /* from j to k */
		    haz = ntrans[jk][i]/nrisk[j][i];
		    htemp = haz/nrisk[j][i];  /* scaled hazard */
		    for (g=0; g <ngrp; g++ )  {
			if (wg[j][g] >0) C[jk][g] -= wg[j][g] * htemp;
		    }

		    if (j != k) {  /* no effect if j==k */
			/* Update U */
			for (g=0; g <ngrp; g++) {
			    if (wg[j][g] > 0) {
				U[j][g] += wg[j][g] * phat[j] * htemp;
				U[k][g] -= wg[j][g] * phat[j]  *htemp;
			    }
			}
		    }
		}
	    }

	    /* update the standard errors */
	    if (tdeath >0 || i==0) {
		for (j=0; j<nstate; j++) {
		    temp=0;
		    for (g=0; g< ngrp; g++) temp += U[j][g]* U[j][g];
		    se1[j] = sqrt(temp);
		}

		for (jk=0; jk<nhaz; jk++) {
		    temp =0;
		    for (g=0; g< ngrp; g++) temp += C[jk][g] * C[jk][g];
		    se2[jk] = sqrt(temp);
		}
	    }
	}   /*end of influence and se compuation */
	
	/* 
	** We can now update phat, p(t-) was needed for the IJ above
	** It is only nstate elements, no need to worry about efficiency
	** If se==0, H will not have been computed, so do it using ntrans
	*/
	for (j=0; j<nstate; j++) hcopy[j] = phat[j]; /* hcopy as a temp */
	for (jk=0; jk < nhaz; jk++){
	    if (ntrans[jk][i] > 0) {
		j = trmat[jk];
		k = trmat[jk + nhaz];
		haz = ntrans[jk][i]/ nrisk[j][i];
		chaz[jk] += haz;
		phat[j] -= hcopy[j]*haz;
		phat[k] += hcopy[j]*haz;
	    }
	}	

	/* Save out the results */
	saveit: j=j+1;  /* dummy line to hang the label on */
	for (j=0; j<nstate; j++) {
	    pstate[j][i] = phat[j];
	    if (sefit >0) {
		stdp[j][i] = se1[j];
		stda[j][i] = se3[j];
	    }
	} 

	for (jk=0; jk< nhaz; jk++) {
	    cumhaz[jk][i] = chaz[jk];
	    if (sefit > 0) stdc[jk][i] = se2[jk];
	}
	
	if (sefit > 1) { /* save the influence matrix */
	    for (j=0; j<nstate; j++) {
		for (g=0; g<ngrp; g++) *usave++ = U[j][g];
	    }
	}
    }
  
    /* return a list */
    UNPROTECT(1);
    return(rlist);
}
