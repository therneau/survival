/*
** A modified version of concordance3 and concordance4: skip the influence
**   and std error portions, for speed.  
*/

#include "survS.h"
#include "survproto.h"
   
SEXP concordance5(SEXP y, SEXP x2, SEXP wt2, SEXP timewt2, 
                      SEXP sortstop) {
    int i, j, k, ii, jj;
    int n, ntree, nevent;
    double *time, *status;
    int xsave;

    /* sum of weights for a node (nwt), sum of weights for the node and
    **  all of its children (twt), then the same again for the subset of
    **  deaths
    */
    double *nwt, *twt;
        
    int ndeath;   /* total number of deaths at this point */    
    int utime;    /* number of unique event times seen so far */
    double dwt, dwt2;   /* sum of weights for deaths and deaths tied on x */
    double wsum[3]; /* the sum of weights that are > current, <, or equal  */
    double adjtimewt;  /* accounts for npair and timewt*/

    SEXP rlist, count2;
    double *count;
    double *wt, *timewt;
    int    *x, *sort2;
    static const char *outnames[]={"count", ""};
      
    n = nrows(y);
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
        
    nwt = (double *) R_alloc(2*ntree, sizeof(double));
    twt = nwt + ntree;
    for (i=0; i< 2*ntree; i++) nwt[i] =0.0;
    
    PROTECT(rlist = mkNamed(VECSXP, outnames));
    count2 = SET_VECTOR_ELT(rlist, 0, allocVector(REALSXP, 5));
    count = REAL(count2); 
    for (i=0; i<5; i++) count[i]=0.0;
    
    utime=0;
    for (i=0; i<n;) {
        ii = sort2[i];  
        if (status[ii]==0) { /* censored, simply add them into the tree */
            addin(nwt, twt, x[ii], wt[ii]);
            i++;
        }
        else {  /* process all tied deaths at this point */
            ndeath=0; dwt=0; 
            dwt2 =0; xsave=x[ii]; 
            adjtimewt = timewt[utime++];

            /* pass 1 */
            for (j=i; j<n && time[sort2[j]]==time[ii]; j++) {
                jj = sort2[j];
                ndeath++; 
                count[3] += wt[jj] * dwt * adjtimewt;  /* update tied on y */
                dwt += wt[jj];   /* sum of wts at this death time */

                if (x[jj] != xsave) {  /* restart the tied.xy counts */
		    dwt2 =0;
                    xsave = x[jj];
                }
                count[4] += wt[jj] * dwt2 * adjtimewt; /* tied on xy */
                dwt2 += wt[jj]; /* sum of tied.xy weights */

                /* Count concordant, discordant, etc. */
                walkup(nwt, twt, x[jj], wsum, ntree);
                for (k=0; k<3; k++) {
                    count[k] += wt[jj]* wsum[k] * adjtimewt;
                }
            }
      
            /* pass 2 */
            for (j=i; j< (i+ndeath); j++) {
                jj = sort2[j];
                addin(nwt, twt, x[jj], wt[jj]); 
            }
            i += ndeath;
        }
    }
    count[3] -= count[4];   /* the tied.xy were counted twice, once as tied.y */
        
    UNPROTECT(1);
    return(rlist);
}
SEXP concordance6(SEXP y, SEXP x2, SEXP wt2, SEXP timewt2, 
                      SEXP sortstart, SEXP sortstop) {
    int i, j, k, ii, jj, i2;
    int n, ntree, nevent;
    double *time1, *time2, *status;
    int xsave; 

    /* sum of weights for a node (nwt), sum of weights for the node and
    **  all of its children (twt), then the same again for the subset of
    **  deaths
    */
    double *nwt, *twt;
        
    int ndeath;   /* total number of deaths at this point */    
    int utime;    /* number of unique event times seen so far */
    double dwt;   /* weighted number of deaths at this point */
    double dwt2;  /* tied on both x and y */
    double wsum[3]; /* the sum of weights that are > current, <, or equal  */
    double adjtimewt;  /* accounts for npair and timewt*/

    SEXP rlist, count2;
    double *count;
    double *wt, *timewt;
    int    *x, *sort2, *sort1;
    static const char *outnames[]={"count", ""};
       
    n = nrows(y);
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
    **  tree holding all subjects. 
    */
    nwt = (double *) R_alloc(2*ntree, sizeof(double));
    twt = nwt + ntree;
    
    for (i=0; i< 4*ntree; i++) nwt[i] =0.0;
    
    PROTECT(rlist = mkNamed(VECSXP, outnames));
    count2 = SET_VECTOR_ELT(rlist, 0, allocVector(REALSXP, 5));
    count = REAL(count2); 
    for (i=0; i<5; i++) count[i]=0.0;
    
    utime=0; i2 =0;  /* i2 tracks the start times */
    for (i=0; i<n;) {
        ii = sort2[i];  
        if (status[ii]==0) { /* censored, simply add them into the tree */
            addin(nwt, twt, x[ii], wt[ii]);
            i++;
        }
        else {  /* a death */
            /* remove any subjects whose start time has been passed */
            for (; i2<n && (time1[sort1[i2]] >= time2[ii]); i2++) {
                jj = sort1[i2];
		addin(nwt, twt, x[jj], -wt[jj]);  /*remove from main tree */
            }

            ndeath=0; dwt=0; 
            dwt2 =0; xsave=x[ii]; 
            adjtimewt = timewt[utime++];

            /* pass 1 */
            for (j=i; j<n && (time2[sort2[j]]==time2[ii]); j++) {
                jj = sort2[j];
                ndeath++; 
                jj = sort2[j];
                count[3] += wt[jj] * dwt;  /* update total tied on y */
                dwt += wt[jj];   /* count of deaths and sum of wts */

                if (x[jj] != xsave) {  /* restart the tied.xy counts */
                    dwt2 =0;
                    xsave = x[jj];
                }
                count[4] += wt[jj] * dwt2 * adjtimewt; /* tied on xy */
                dwt2 += wt[jj]; /* sum of tied.xy weights */

                /* Count concordant, discordant, etc. */
                walkup(nwt, twt, x[jj], wsum, ntree);
                for (k=0; k<3; k++) {
                    count[k] += wt[jj]* wsum[k] * adjtimewt;
		}
            }
            /* pass 2 */
            for (j=i; j< (i+ndeath); j++) {
                jj = sort2[j];
                addin(nwt, twt, x[jj], wt[jj]); 
            }

            i += ndeath;

        }
    }

    count[3] -= count[4]; /* tied.y was double counted a tied.xy */
        
    UNPROTECT(1);
    return(rlist);
}
