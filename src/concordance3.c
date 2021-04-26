/* Automatically generated from the noweb directory */
#include "survS.h"

void walkup(double *nwt, double* twt, int index, double sums[3], int ntree) {
    int i, j, parent;

    for (i=0; i<3; i++) sums[i] = 0.0;
    sums[2] = nwt[index];   /* tied on x */
    
    j = 2*index +2;  /* right child */
    if (j < ntree) sums[0] += twt[j];
    if (j <=ntree) sums[1]+= twt[j-1]; /*left child */

    while(index > 0) { /* for as long as I have a parent... */
        parent = (index-1)/2;
        if (index%2 == 1) sums[0] += twt[parent] - twt[index]; /* left child */
        else sums[1] += twt[parent] - twt[index]; /* I am a right child */
        index = parent;
    }
}

void addin(double *nwt, double *twt, int index, double wt) {
    nwt[index] += wt;
    while (index >0) {
        twt[index] += wt;
        index = (index-1)/2;
    }
    twt[0] += wt;
}
    
SEXP concordance3(SEXP y, SEXP x2, SEXP wt2, SEXP timewt2, 
                      SEXP sortstop, SEXP doresid2) {
    int i, j, k, ii, jj, kk, j2;
    int n, ntree, nevent;
    double *time, *status;
    int xsave;

    /* sum of weights for a node (nwt), sum of weights for the node and
    **  all of its children (twt), then the same again for the subset of
    **  deaths
    */
    double *nwt, *twt, *dnwt, *dtwt;
    double z2;  /* sum of z^2 values */    
        
    int ndeath;   /* total number of deaths at this point */    
    int utime;    /* number of unique event times seen so far */
    double dwt, dwt2;   /* sum of weights for deaths and deaths tied on x */
    double wsum[3]; /* the sum of weights that are > current, <, or equal  */
    double temp, adjtimewt;  /* the second accounts for npair and timewt*/

    SEXP rlist, count2, imat2, resid2;
    double *count, *imat[5], *resid[4];
    double *wt, *timewt;
    int    *x, *sort2;
    int doresid;
    static const char *outnames1[]={"count", "influence", "resid", ""},
                      *outnames2[]={"count", "influence", ""};
      
    n = nrows(y);
    doresid = asLogical(doresid2);
    x = INTEGER(x2);
    wt = REAL(wt2);
    timewt = REAL(timewt2);
    sort2 = INTEGER(sortstop);
    time = REAL(y);
    status = time + n;
   
    /* if there are tied predictors, the total size of the tree will be < n */
    ntree =0; nevent =0;
    for (i=0; i<n; i++) {
        if (x[i] >= ntree) ntree = x[i] +1;  
        nevent += status[i];
    }
        
    nwt = (double *) R_alloc(4*ntree, sizeof(double));
    twt = nwt + ntree;
    dnwt = twt + ntree;
    dtwt = dnwt + ntree;
    
    for (i=0; i< 4*ntree; i++) nwt[i] =0.0;
    
    if (doresid) PROTECT(rlist = mkNamed(VECSXP, outnames1));
    else  PROTECT(rlist = mkNamed(VECSXP, outnames2));
    count2 = SET_VECTOR_ELT(rlist, 0, allocVector(REALSXP, 6));
    count = REAL(count2); 
    for (i=0; i<6; i++) count[i]=0.0;
    imat2 = SET_VECTOR_ELT(rlist, 1, allocMatrix(REALSXP, n, 5));
    for (i=0; i<5; i++) {
        imat[i] = REAL(imat2) + i*n;
        for (j=0; j<n; j++) imat[i][j] =0;
    }
    if (doresid==1) {
        resid2 = SET_VECTOR_ELT(rlist, 2, allocMatrix(REALSXP, nevent, 4));
        for (i=0; i<4; i++) resid[i] = REAL(resid2) + i*nevent;
        }
    
    z2 =0; utime=0;
    for (i=0; i<n;) {
        ii = sort2[i];  
        if (status[ii]==0) { /* censored, simply add them into the tree */
            /* Initialize the influence */
            walkup(dnwt, dtwt, x[ii], wsum, ntree);
            imat[0][ii] -= wsum[1];
            imat[1][ii] -= wsum[0];
            imat[2][ii] -= wsum[2];
            
            /* Cox variance */
            walkup(nwt, twt, x[ii], wsum, ntree);
            z2 += wt[ii]*(wsum[0]*(wt[ii] + 2*(wsum[1] + wsum[2])) +
                          wsum[1]*(wt[ii] + 2*(wsum[0] + wsum[2])) +
                          (wsum[0]-wsum[1])*(wsum[0]-wsum[1]));
            /* add them to the tree */
            addin(nwt, twt, x[ii], wt[ii]);
            i++;
        }
        else {  /* process all tied deaths at this point */
            ndeath=0; dwt=0; 
            dwt2 =0; xsave=x[ii]; j2= i;
            adjtimewt = timewt[utime++];

            /* pass 1 */
            for (j=i; j<n && time[sort2[j]]==time[ii]; j++) {
                jj = sort2[j];
                ndeath++; 
                count[3] += wt[jj] * dwt * adjtimewt;  /* update total tied on y */
                dwt += wt[jj];   /* sum of wts at this death time */

                if (x[jj] != xsave) {  /* restart the tied.xy counts */
                    if (wt[sort2[j2]] < dwt2) { /* more than 1 tied */
                        for (; j2<j; j2++) {
                            /* update influence for this subgroup of x */
                            kk = sort2[j2];
                            imat[4][kk] += (dwt2- wt[kk]) * adjtimewt;
                            imat[3][kk] -= (dwt2- wt[kk]) * adjtimewt;
                        }
                    } else j2 = j;
                    dwt2 =0;
                    xsave = x[jj];
                }
                count[4] += wt[jj] * dwt2 * adjtimewt; /* tied on xy */
                dwt2 += wt[jj]; /* sum of tied.xy weights */

                /* Count concordant, discordant, etc. */
                walkup(nwt, twt, x[jj], wsum, ntree);
                for (k=0; k<3; k++) {
                    count[k] += wt[jj]* wsum[k] * adjtimewt;
                    imat[k][jj] += wsum[k]*adjtimewt;
                }

                /* add to the event tree */
                addin(dnwt, dtwt, x[jj], adjtimewt*wt[jj]);  /* weighted deaths */

                /* first part of residuals */
                if (doresid) {
                    nevent--;
                    resid[0][nevent] = (wsum[0] - wsum[1])/twt[0]; /* -1 to 1 */
                    resid[1][nevent] = twt[0] * adjtimewt;
                    resid[2][nevent] = wt[jj];
                }
            }
            /* finish the tied.xy influence */
            if (wt[sort2[j2]] < dwt2) { /* more than 1 tied */
                for (; j2<j; j2++) {
                    /* update influence for this subgroup of x */
                    kk = sort2[j2];
                    imat[4][kk] += (dwt2- wt[kk]) * adjtimewt;
                    imat[3][kk] -= (dwt2- wt[kk]) * adjtimewt;
                }
            }
      
            /* pass 2 */
            for (j=i; j< (i+ndeath); j++) {
                jj = sort2[j];
                /* Update influence */
                walkup(dnwt, dtwt, x[jj], wsum, ntree);
                imat[0][jj] -= wsum[1];
                imat[1][jj] -= wsum[0];
                imat[2][jj] -= wsum[2];  /* tied.x */
                imat[3][jj] += (dwt- wt[jj])* adjtimewt;
     
                /* increment Cox var and add obs into the tree */
                walkup(nwt, twt, x[jj], wsum, ntree);
                z2 += wt[jj]*(wsum[0]*(wt[jj] + 2*(wsum[1] + wsum[2])) +
                              wsum[1]*(wt[jj] + 2*(wsum[0] + wsum[2])) +
                              (wsum[0]-wsum[1])*(wsum[0]-wsum[1]));

                addin(nwt, twt, x[jj], wt[jj]); 
            }
            count[5] += dwt * adjtimewt* z2/twt[0]; /* weighted var in risk set*/
            i += ndeath;

            if (doresid) { /*Add the last part of the residuals */
                temp = twt[0]*twt[0]*twt[0];
                for (j=0; j<ndeath; j++)
                    resid[3][nevent+j] = z2/temp;
            }
        }
    }

    /* 
    ** Now finish off the influence for each observation 
    **  Since times flip (looking backwards) the wsum contributions flip too
    */
    for (i=0; i<n; i++) {
        ii = sort2[i];
        walkup(dnwt, dtwt, x[ii], wsum, ntree);
        imat[0][ii] += wsum[1];
        imat[1][ii] += wsum[0];
        imat[2][ii] += wsum[2];
    }
    count[3] -= count[4];   /* the tied.xy were counted twice, once as tied.y */
        
    UNPROTECT(1);
    return(rlist);
}
    SEXP concordance4(SEXP y, SEXP x2, SEXP wt2, SEXP timewt2, 
                      SEXP sortstart, SEXP sortstop, SEXP doresid2) {
    int i, j, k, ii, jj, kk, i2, j2;
    int n, ntree, nevent;
    double *time1, *time2, *status;
    int xsave; 

    /* sum of weights for a node (nwt), sum of weights for the node and
    **  all of its children (twt), then the same again for the subset of
    **  deaths
    */
    double *nwt, *twt, *dnwt, *dtwt;
    double z2;  /* sum of z^2 values */    
        
    int ndeath;   /* total number of deaths at this point */    
    int utime;    /* number of unique event times seen so far */
    double dwt;   /* weighted number of deaths at this point */
    double dwt2;  /* tied on both x and y */
    double wsum[3]; /* the sum of weights that are > current, <, or equal  */
    double temp, adjtimewt;  /* the second accounts for npair and timewt*/

    SEXP rlist, count2, imat2, resid2;
    double *count, *imat[5], *resid[4];
    double *wt, *timewt;
    int    *x, *sort2, *sort1;
    int doresid;
    static const char *outnames1[]={"count", "influence", "resid", ""},
                      *outnames2[]={"count", "influence", ""};
      
    n = nrows(y);
    doresid = asLogical(doresid2);
    x = INTEGER(x2);
    wt = REAL(wt2);
    timewt = REAL(timewt2);
    sort2 = INTEGER(sortstop);
    sort1 = INTEGER(sortstart);
    time1 = REAL(y);
    time2 = time1 + n;
    status = time2 + n;
   
    /* if there are tied predictors, the total size of the tree will be < n */
    ntree =0; nevent =0;
    for (i=0; i<n; i++) {
        if (x[i] >= ntree) ntree = x[i] +1;  
        nevent += status[i];
    }
        
    /*
    ** nwt and twt are the node weight and total =node + all children for the
    **  tree holding all subjects.  dnwt and dtwt are the same for the tree
    **  holding all the events
    */
    nwt = (double *) R_alloc(4*ntree, sizeof(double));
    twt = nwt + ntree;
    dnwt = twt + ntree;
    dtwt = dnwt + ntree;
    
    for (i=0; i< 4*ntree; i++) nwt[i] =0.0;
    
    if (doresid) PROTECT(rlist = mkNamed(VECSXP, outnames1));
    else  PROTECT(rlist = mkNamed(VECSXP, outnames2));
    count2 = SET_VECTOR_ELT(rlist, 0, allocVector(REALSXP, 6));
    count = REAL(count2); 
    for (i=0; i<6; i++) count[i]=0.0;
    imat2 = SET_VECTOR_ELT(rlist, 1, allocMatrix(REALSXP, n, 5));
    for (i=0; i<5; i++) {
        imat[i] = REAL(imat2) + i*n;
        for (j=0; j<n; j++) imat[i][j] =0;
    }
    if (doresid==1) {
        resid2 = SET_VECTOR_ELT(rlist, 2, allocMatrix(REALSXP, nevent, 4));
        for (i=0; i<4; i++) resid[i] = REAL(resid2) + i*nevent;
        }
    
    z2 =0; utime=0; i2 =0;  /* i2 tracks the start times */
    for (i=0; i<n;) {
        ii = sort2[i];  
        if (status[ii]==0) { /* censored, simply add them into the tree */
            /* Initialize the influence */
            walkup(dnwt, dtwt, x[ii], wsum, ntree);
            imat[0][ii] -= wsum[1];
            imat[1][ii] -= wsum[0];
            imat[2][ii] -= wsum[2];
            
            /* Cox variance */
            walkup(nwt, twt, x[ii], wsum, ntree);
            z2 += wt[ii]*(wsum[0]*(wt[ii] + 2*(wsum[1] + wsum[2])) +
                          wsum[1]*(wt[ii] + 2*(wsum[0] + wsum[2])) +
                          (wsum[0]-wsum[1])*(wsum[0]-wsum[1]));
            /* add them to the tree */
            addin(nwt, twt, x[ii], wt[ii]);
            i++;
        }
        else {  /* a death */
            /* remove any subjects whose start time has been passed */
            for (; i2<n && (time1[sort1[i2]] >= time2[ii]); i2++) {
                jj = sort1[i2];
                /* influence */
                walkup(dnwt, dtwt, x[jj], wsum, ntree);
                imat[0][jj] += wsum[1];
                imat[1][jj] += wsum[0];
                imat[2][jj] += wsum[2];

                addin(nwt, twt, x[jj], -wt[jj]);  /*remove from main tree */

                /* Cox variance */
                walkup(nwt, twt, x[jj], wsum, ntree);
                z2 -= wt[jj]*(wsum[0]*(wt[jj] + 2*(wsum[1] + wsum[2])) +
                              wsum[1]*(wt[jj] + 2*(wsum[0] + wsum[2])) +
                              (wsum[0]-wsum[1])*(wsum[0]-wsum[1]));
            }

            ndeath=0; dwt=0; 
            dwt2 =0; xsave=x[ii]; j2= i;
            adjtimewt = timewt[utime++];

            /* pass 1 */
            for (j=i; j<n && (time2[sort2[j]]==time2[ii]); j++) {
                jj = sort2[j];
                ndeath++; 
                jj = sort2[j];
                count[3] += wt[jj] * dwt;  /* update total tied on y */
                dwt += wt[jj];   /* count of deaths and sum of wts */

                if (x[jj] != xsave) {  /* restart the tied.xy counts */
                    if (wt[sort2[j2]] < dwt2) { /* more than 1 tied */
                        for (; j2<j; j2++) {
                            /* update influence for this subgroup of x */
                            kk = sort2[j2];
                            imat[4][kk] += (dwt2- wt[kk]) * adjtimewt;
                            imat[3][kk] -= (dwt2- wt[kk]) * adjtimewt;
                        }
                    } else j2 = j;
                    dwt2 =0;
                    xsave = x[jj];
                }
                count[4] += wt[jj] * dwt2 * adjtimewt; /* tied on xy */
                dwt2 += wt[jj]; /* sum of tied.xy weights */

                /* Count concordant, discordant, etc. */
                walkup(nwt, twt, x[jj], wsum, ntree);
                for (k=0; k<3; k++) {
                    count[k] += wt[jj]* wsum[k] * adjtimewt;
                    imat[k][jj] += wsum[k]*adjtimewt;
                }

                /* add to the event tree */
                addin(dnwt, dtwt, x[jj], adjtimewt*wt[jj]);  /* weighted deaths */

                /* first part of residuals */
                if (doresid) {
                    nevent--;
                    resid[0][nevent] = (wsum[0] - wsum[1])/twt[0]; /* -1 to 1 */
                    resid[1][nevent] = twt[0] * adjtimewt;
                    resid[2][nevent] = wt[jj];
                }
            }
            /* finish the tied.xy influence */
            if (wt[sort2[j2]] < dwt2) { /* more than 1 tied */
                for (; j2<j; j2++) {
                    /* update influence for this subgroup of x */
                    kk = sort2[j2];
                    imat[4][kk] += (dwt2- wt[kk]) * adjtimewt;
                    imat[3][kk] -= (dwt2- wt[kk]) * adjtimewt;
                }
            }

            /* pass 3 */
            for (j=i; j< (i+ndeath); j++) {
                jj = sort2[j];
                /* Update influence */
                walkup(dnwt, dtwt, x[jj], wsum, ntree);
                imat[0][jj] -= wsum[1];
                imat[1][jj] -= wsum[0];
                imat[2][jj] -= wsum[2];  /* tied.x */
                imat[3][jj] += (dwt- wt[jj])* adjtimewt;

                /* increment Cox var and add obs into the tree */
                walkup(nwt, twt, x[jj], wsum, ntree);
                z2 += wt[jj]*(wsum[0]*(wt[jj] + 2*(wsum[1] + wsum[2])) +
                              wsum[1]*(wt[jj] + 2*(wsum[0] + wsum[2])) +
                              (wsum[0]-wsum[1])*(wsum[0]-wsum[1]));

                addin(nwt, twt, x[jj], wt[jj]); 
            }
            count[5] += dwt * adjtimewt* z2/twt[0]; /* weighted var in risk set*/
            i += ndeath;

            if (doresid) { /*Add the last part of the residuals */
                temp = twt[0]*twt[0]*twt[0];
                for (j=0; j<ndeath; j++)
                    resid[3][nevent+j] = z2/temp;
            }
        }
    }

    /* 
    ** Now finish off the influence for those not yet removed
    **  Since times flip (looking backwards) the wsum contributions flip too
    */
    for (; i2<n; i2++) {
        ii = sort1[i2];
        walkup(dnwt, dtwt, x[ii], wsum, ntree);
        imat[0][ii] += wsum[1];
        imat[1][ii] += wsum[0];
        imat[2][ii] += wsum[2];
    }
    count[3] -= count[4]; /* tied.y was double counted a tied.xy */
        
    UNPROTECT(1);
    return(rlist);
}
