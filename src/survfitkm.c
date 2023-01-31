/* Automatically generated from the noweb directory */
#include <math.h>
#include "survS.h"
#include "survproto.h"

SEXP survfitkm(SEXP y2, SEXP weight2,  SEXP sort12, SEXP sort22, 
               SEXP type2, SEXP id2, SEXP nid2,   SEXP position2,
               SEXP influence2) {
              
    int i, i1, i2, j, k, person1, person2;
    int nused, nid, type, influence;
    int ny, ntime;
    double *tstart=0, *stime, *status, *wt;
    double v1, v2, dtemp, haz;
    double temp, dtemp2, dtemp3, frac, btemp;
    int *sort1=0, *sort2, *id=0;
    static const char *outnames[]={"time", "n", "estimate", "std.err",
                                     "influence1", "influence2", ""};
    SEXP rlist;
    double *gwt=0, *inf1=0, *inf2=0;  /* =0 to silence -Wall */
    int *gcount=0;
    int n1, n2, n3, n4;
    int *position=0, hasid;
    double wt1, wt2, wt3, wt4;
                      
    /* output variables */
    double  *n[8],  *dtime,
            *kvec, *nvec, *std[2], *imat1=0, *imat2=0; /* =0 to silence -Wall*/
    double km, nelson;  /* current estimates */

    /* map the input data */
    ny = ncols(y2);     /* 2= ordinary survival 3= start,stop data */
    nused = nrows(y2);
    if (ny==3) { 
        tstart = REAL(y2);
        stime = tstart + nused;
        sort1 = INTEGER(sort12);
    }
    else stime = REAL(y2);
    status= stime +nused;
    wt = REAL(weight2);
    sort2 = INTEGER(sort22);
    nused = LENGTH(sort22);
                   
    type = asInteger(type2);
    nid = asInteger(nid2);
    if (LENGTH(position2) > 0) {
        hasid =1;
        position = INTEGER(position2);
    } else hasid=0;
    influence = asInteger(influence2);

    /* nused was used for two things just above.  The first was the length of
       the input data y, only needed for a moment to set up tstart, stime, and
       status.  The second is the number of these observations we will actually
       use, which is the length of sort2.  This routine can be called multiple
       times with sort1/sort2 pointing to different subsets of the data while
       y, wt, id and position can remain unchanged
    */

    if (length(id2)==0) nid =0;  /* no robust variance */
    else id = INTEGER(id2);

    /* pass 1, get the number of unique times, needed for memory allocation 
      Number of xval groups (unique id values) has been supplied 
      Data is sorted by time
    */
    ntime =1; 
    temp = stime[sort2[0]];
    for (i=1; i<nused; i++) {
        i2 = sort2[i];
        if (stime[i2] != temp) {
            ntime++;
            temp = stime[i2];
        }
    }        
    /* Allocate memory for the output 
        n has 6 columns for number at risk, events, censor, then the 
        3 weighted versions of the same, then optionally two more for
        number added to the risk set (when ny=3)
    */
    PROTECT(rlist = mkNamed(VECSXP, outnames));
    
    dtime  = REAL(SET_VECTOR_ELT(rlist, 0, allocVector(REALSXP, ntime)));
    if (ny==2) j=6;  else j=8;
    n[0]  = REAL(SET_VECTOR_ELT(rlist, 1, allocMatrix(REALSXP, ntime, j)));
    for (i=1; i<j; i++) n[i] = n[0] + i*ntime;

    kvec  = REAL(SET_VECTOR_ELT(rlist, 2, allocMatrix(REALSXP, ntime, 2)));
    nvec  = kvec + ntime;  /* Nelson-Aalen estimate */
    std[0] = REAL(SET_VECTOR_ELT(rlist, 3, allocMatrix(REALSXP, ntime,2)));
    std[1] = std[0] + ntime;
  
    if (nid >0 ) { /* robust variance */
        gcount = (int *) R_alloc(nid, sizeof(int));
        if (type <3) {  /* working vectors for the influence */
            gwt  = (double *) R_alloc(3*nid, sizeof(double)); 
            inf1 = gwt + nid;
            inf2 = inf1 + nid; 
            for (i=0; i< nid; i++) {
                gwt[i] =0.0;
                gcount[i] = 0;
                inf1[i] =0;
                inf2[i] =0;
            }
        }
        else {
            gwt = (double *) R_alloc(2*nid, sizeof(double));
            inf2 = gwt + nid;
            for (i=0; i< nid; i++) {
                gwt[i] =0.0;
                gcount[i] = 0;
                inf2[i] =0;
            }
        }

        /* these are not accumulated, so do not need to be zeroed */
        if (type <3) {
            if (influence==1 || influence ==3) 
                imat1 = REAL(SET_VECTOR_ELT(rlist, 4,
                                     allocMatrix(REALSXP, nid, ntime)));
            if (influence==2 || influence==3) 
                imat2 =  REAL(SET_VECTOR_ELT(rlist, 5,
                           allocMatrix(REALSXP, nid, ntime))); 
        }                
        else  if (influence !=0) 
                imat2 = REAL(SET_VECTOR_ELT(rlist, 5,
                               allocMatrix(REALSXP, nid, ntime))); 
    }

    R_CheckUserInterrupt();  /*check for control-C */
    /*
    ** person1, person2 track through sort1 and sort2, respectively
    **  likewise with i1 and i2
    */
    person1 = nused-1;  person2 = nused-1;
    n1=0; wt1=0;
    for (k=ntime-1; k>=0; k--) {
        dtime[k] = stime[sort2[person2]];  /* current time point */
        n2=0; n3=0; wt2=0; wt3=0; 
        for (; person2 >=0; person2--) {
            i2= sort2[person2];
            if (stime[i2] != dtime[k]) break;

            n1++;             /* number at risk */
            wt1 += wt[i2];    /* weighted number at risk */
            if (status[i2] ==1) {
                n2++;            /* events */
                wt2 += wt[i2];
            } else if (hasid==0 || (position[i2]& 2)) {
                /* if there are no repeated obs for a subject (hasid=0)
                **  or this is the last of a string (a,b](b,c](c,d].. for
                **  a subject (position[i2]=2 or 3), then it is a 'real' censor
                */
                n3++;
                wt3 += wt[i2];
            }
        }
        
        if (ny==3) { /* remove any with start time >=dtime*/
            n4 =0; wt4 =0;
            for (; person1 >=0; person1--) {
                i1 = sort1[person1];
                if (tstart[i1] < dtime[k]) break;
                n1--;
                wt1 -= wt[i1];
                if (hasid==0 || (position[i1] & 1)) {
                    /* if there are no repeated id (hasid=0) or this is the
                    ** first of a string of (a,b](b,c](c,d] for a subject, then
                    ** this is a 'real' entry */
                    n4++;
                    wt4 += wt[i1];
                }
                n[6][k] = n4;
                n[7][k] = wt4;
            }
        }

        n[0][k] = n1;  n[1][k]=n2;  n[2][k]=n3;
        n[3][k] = wt1; n[4][k]=wt2; n[5][k]=wt3; 
    }

    if (ny ==3) {   /* fill in number entered for the initial interval */
        n4=0; wt4=0;
        for (; person1>=0; person1--) {
            i1 = sort1[person1];
            if (hasid==0 || (position[i1] & 1)) {
                n4++;
                wt4 += wt[i1];
            }
        }
        n[6][0] = n4;    
        n[7][0] = wt4;
    }
    R_CheckUserInterrupt();  /*check for control-C */
    nelson =0.0; km=1.0; 
    v1=0; v2=0;
    if (nid==0) {  /* simple variance */
       if (type==1 || type==3) {  /* Nelson-Aalen hazard */
            for (i=0; i<ntime; i++) {
                if (n[1][i]>0 && n[4][i]>0) {  /* at least one event with wt>0*/
                    nelson += n[4][i]/n[3][i];
                    v2 += n[4][i]/(n[3][i]*n[3][i]);
                    }

                nvec[i] = nelson;
                std[0][i] = sqrt(v2);
                std[1][i] = sqrt(v2);
                }
       } else {              /* Fleming hazard */
            for (i=0; i<ntime; i++) {
                for (j=0; j<n[1][i]; j++) {
                    dtemp = n[3][i] - j*n[4][i]/n[1][i];
                    nelson += n[4][i] /(n[1][i]* dtemp);
                    v2 += n[4][i]/(n[1][i]*dtemp*dtemp);
                }
                kvec[i] = exp(-nelson);
                nvec[i] = nelson;
                std[0][i] = sqrt(v2);
                std[1][i] = sqrt(v2);
            }
        }

        if (type < 3) {  /* KM survival */
            for (i=0; i<ntime; i++) {
                if (n[1][i]>0 && n[4][i]>0) {  /* at least one event */
                    km *= (n[3][i]-n[4][i])/n[3][i];
                    v1 += n[4][i]/(n[3][i] * (n[3][i] - n[4][i])); /* Greenwood */
                    }
                kvec[i] = km;
                std[0][i] = sqrt(v1);
            }
        } else {  /* exp survival */
            for (i=0; i< ntime; i++) {
                kvec[i] = exp(-nvec[i]);
                std[0][i] = std[1][i];
            }
        }
    }

    else { /* infinitesimal jackknife variance */
        v1=0; v2 =0; km=1; nelson =0;
        person2=0; 
        if (ny==3) {
            person1 =0;
        } else {
            /* at the start, everyone is at risk */
            for (i=0; i< nused; i++) {
                i2 = id[i];
                gcount[i2]++;
                gwt[i2] += wt[i];
            }
        }
            
        if (type==1) {
            for (i=0; i< ntime; i++) {
                if (ny==3) {
                    /* add in new subjects */
                    for (; person1 < nused; person1++) {
                        /* add in those whose start time is < dtime */
                        i1 = sort1[person1];
                        if (tstart[i1] >= dtime[i]) break;  
                        gcount[id[i1]]++;
                        gwt[id[i1]] += wt[i1];
                    }
                }
         
                if (n[1][i] > 0 && n[4][i]>0) { /* need to update the sums */
                    haz = n[4][i]/n[3][i];
                    for (k=0; k< nid; k++) {
                        inf1[k] = inf1[k] *(1.0 -haz) + gwt[k]*km*haz/n[3][i];
                        inf2[k] -= gwt[k] * haz/n[3][i];
                    }
                    for (; person2<nused; person2++) { 
                        i2 = sort2[person2];
                        if (stime[i2] > dtime[i]) break;   /* those at this time */
                        if (status[i2]==1) {
                            inf1[id[i2]] -= km* wt[i2]/n[3][i];
                            inf2[id[i2]] += wt[i2]/n[3][i];
                        }
                        gcount[id[i2]] --;
                        if (gcount[id[i2]] ==0) gwt[id[i2]] = 0.0;
                        else gwt[id[i2]] -= wt[i2];
                    }
                    km *= (1-haz);
                    nelson += haz;
                   
                    v1=0; v2=0;
                    for (k=0; k<nid; k++) {
                        v1 += inf1[k]*inf1[k];
                        v2 += inf2[k]*inf2[k];
                    }
                } else {  /* only need to udpate weights */
                    for (; person2<nused; person2++) { 
                        i2 = sort2[person2];
                        if (stime[i2] > dtime[i]) break;
                        gcount[id[i2]] --;
                        if (gcount[id[i2]] ==0) gwt[id[i2]] = 0.0;
                        else gwt[id[i2]] -= wt[i2];
                    }
                }
         
                kvec[i] = km;
                nvec[i] = nelson;
                std[0][i] = sqrt(v1);
                std[1][i] = sqrt(v2);
                if (influence==1 || influence ==3) 
                    for (k=0; k<nid; k++) *imat1++ = inf1[k];
                if (influence==2 || influence ==3)
                    for (k=0; k<nid; k++) *imat2++ = inf2[k];
            }
        }
        else if (type==2) {  /* KM survival, Fleming-Harrington hazard */
            for (i=0; i< ntime; i++) {
                if (ny==3) {
                    /* add in new subjects */
                    for (; person1<nused; person1++) {
                        i1 = sort1[person1];
                        if (tstart[i1] >= dtime[i]) break;
                        gcount[id[i1]]++;
                        gwt[id[i1]] += wt[i1];
                    }
                }
         
                if (n[1][i] > 0 && n[4][i] >0) { /* need to update the sums */
                    dtemp =0;  /* the working denominator */
                    dtemp2=0;  /* sum of squares */
                    dtemp3=0;
                    temp = n[3][i] - n[4][i];  /* sum of weights for the non-deaths */
                    for (k=n[1][i]; k>0; k--) {
                        frac = k/n[1][i];
                        btemp = 1/(temp + frac*n[4][i]);  /* "b" in the math */
                        dtemp += btemp;
                        dtemp2 += btemp*btemp*frac;
                        dtemp3 += btemp*btemp;    /* non-death deriv */
                    }

                    dtemp /=  n[1][i];        /* average denominator */
                    if (n[4][i] != n[1][i]) { /* case weights */
                        dtemp2 *= n[4][i]/ n[1][i];
                        dtemp3 *= n[4][i]/ n[1][i];
                    }
                    nelson += n[4][i]*dtemp;

                    haz = n[4][i]/n[3][i];
                    for (k=0; k< nid; k++) {
                        inf1[k] = inf1[k] *(1.0 -haz) + gwt[k]*km*haz/n[3][i];
                        if (gcount[k]>0) inf2[k] -= gwt[k] * dtemp3;
                    }
                    for (; person2<nused; person2++) {                
                        /* catch the endpoints up to this event time */
                        i2 = sort2[person2];
                        if (stime[i2] > dtime[i]) break;
                        if (status[i2]==1) {
                            inf1[id[i2]] -= km* wt[i2]/n[3][i];
                            inf2[id[i2]] += wt[i2] *(dtemp + dtemp3 - dtemp2);
                         }
                        gcount[id[i2]] --;
                        if (gcount[id[i2]] ==0) gwt[id[i2]] = 0.0;
                        else gwt[id[i2]] -= wt[i2];
                    }
                    km *= (1-haz);
                    
                    v1=0; v2=0;
                    for (k=0; k<nid; k++) {
                        v1 += inf1[k]*inf1[k];
                        v2 += inf2[k]*inf2[k];
                    }
                } else {  /* only need to udpate weights */
                    for (; person2<nused; person2++) { 
                        i2 = sort2[person2];
                        if (stime[i2] > dtime[i]) break;
                        gcount[id[i2]] --;
                        if (gcount[id[i2]] ==0) gwt[id[i2]] = 0.0;
                        else gwt[id[i2]] -= wt[i2];
                    }
                }
         
                kvec[i] = km;
                nvec[i] = nelson;
                std[0][i] = sqrt(v1);
                std[1][i] = sqrt(v2);
                if (influence==1 || influence ==3) 
                    for (k=0; k<nid; k++) *imat1++ = inf1[k];
                if (influence==2 || influence ==3)
                    for (k=0; k<nid; k++) *imat2++ = inf2[k];
            }
        }

        else if (type==3) {  /* exp() survival, NA hazard */
            for (i=0; i< ntime; i++) {
                if (ny==3) {
                    /* add in new subjects */
                    for (; person1 < nused; person1++) {
                        i1 = sort1[person1];
                        if (tstart[i1] >= dtime[i]) break;

                        gcount[id[i1]]++;
                        gwt[id[i1]] += wt[i1];
                    }
                }
                 
                if (n[1][i] > 0 && n[4][i]>0) { /* need to update the sums */
                    haz = n[4][i]/n[3][i];
                    for (k=0; k< nid; k++) {
                        inf2[k] -= gwt[k] * haz/n[3][i];
                    }
                    for (; person2<nused; person2++) { 
                        /* catch the endpoints up to this event time */
                        i2 = sort2[person2];
                        if (stime[i2] > dtime[i]) break;
                        if (status[i2]==1) {
                             inf2[id[i2]] += wt[i2]/n[3][i];
                        }
                        gcount[id[i2]] --;
                        if (gcount[id[i2]] ==0) gwt[id[i2]] = 0.0;
                        else gwt[id[i2]] -= wt[i2];
                    }
                    nelson += haz;
                   
                    v2=0;
                    for (k=0; k<nid; k++) {
                        v2 += inf2[k]*inf2[k];
                    }
                } else {  /* only need to udpate weights */
                    for (; person2<nused; person2++) { 
                        i2 = sort2[person2];
                        if (stime[i2] > dtime[i]) break;
                        gcount[id[i2]] --;
                        if (gcount[id[i2]] ==0) gwt[id[i2]] = 0.0;
                        else gwt[id[i2]] -= wt[i2];
                    }
                }
         
                kvec[i] = exp(-nelson);
                nvec[i] = nelson;
                std[1][i] = sqrt(v2);
                std[0][i] = sqrt(v2);
                
                if (influence>0)
                    for (k=0; k<nid; k++) *imat2++ = inf2[k];
            }

        } else {  /* exp() survival,  Fleming-Harrington hazard */
            for (i=0; i< ntime; i++) {
                if (ny==3) {
                    /* add in new subjects */
                    for (; person1 < nused; person1++) {
                        i1 = sort1[person1];
                        if (tstart[i1] >= dtime[i]) break;
                        gcount[id[i1]]++;
                        gwt[id[i1]] += wt[i1];
                    }
                }
                if (n[1][i] > 0 && n[4][i] >0) { /* need to update the sums */
                        dtemp =0;  /* the working denominator */
                    dtemp2=0;  /* sum of squares */
                    dtemp3=0;
                    temp = n[3][i] - n[4][i];  /* sum of weights for the non-deaths */
                    for (k=n[1][i]; k>0; k--) {
                        frac = k/n[1][i];  
                        btemp = 1/(temp + frac*n[4][i]);  /* "b" in the math */
                        dtemp += btemp;
                        dtemp2 += btemp*btemp*frac;
                        dtemp3 += btemp*btemp;    /* non-death deriv */
                    } 
                    
                    dtemp /=  n[1][i];        /* average denominator */
                    if (n[4][i] != n[1][i]) { /* case weights */
                        dtemp2 *= n[4][i]/ n[1][i];
                        dtemp3 *= n[4][i]/ n[1][i];
                    }
                    nelson += n[4][i]*dtemp;

                    for (k=0; k< nid; k++) {
                        if (gcount[k]>0) inf2[k] -= gwt[k] * dtemp3;
                    }
                    for (; person2<nused; person2++) { 
                         i2 = sort2[person2];
                        if (stime[i2] > dtime[i]) break;
                        if (status[i2]==1) {
                            inf2[id[i2]] += wt[i2] *(dtemp + dtemp3 - dtemp2);
                        }
                        gcount[id[i2]] --;
                        if (gcount[id[i2]] ==0) gwt[id[i2]] = 0.0;
                        else gwt[id[i2]] -= wt[i2];
                    }
            
                    v2=0;
                    for (k=0; k<nid; k++) v2 += inf2[k]*inf2[k];
                }
                else { /* only need to update weights */
                    for (; person2<nused; person2++) { 
                        i2 = sort2[person2];
                        if (stime[i2] > dtime[i]) break;
                        gcount[id[i2]] --;
                        if (gcount[id[i2]] ==0) gwt[id[i2]] = 0.0;
                        else gwt[id[i2]] -= wt[i2];
                    }
                }
                
                kvec[i] = exp(-nelson);
                nvec[i] = nelson;
                std[1][i] = sqrt(v2);
                std[0][i] = sqrt(v2);
                
                if (influence>0)
                    for (k=0; k<nid; k++) *imat2++ = inf2[k];
            }
        }
    }
    
    UNPROTECT(1);
    return(rlist);
}
