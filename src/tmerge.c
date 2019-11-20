/*
** Fill in a new time-dependent variable
** 
**  id: subjec id in the baseline datata set (integer)
**  time1, time2: time intervals in the baseline data, only the second is needed
**  nid, ntime: the id and time point for the new covariate
**  x:  vaule of the new covariate
**  indx:  starting point for matching in the baseline
**
**  the baseline need not be in sorted id order, but all time points
**   for a given id are together, and in order.
*/
#include "survS.h"
#include "survproto.h"

SEXP tmerge(SEXP id2,  SEXP time2x, SEXP newx2,
            SEXP nid2, SEXP ntime2, SEXP x2,  SEXP indx2) {
    int i, k;
    int n1, n2;

    int *id, *nid;
    double *time2, 
	   *ntime, *x;
    int *indx;
    double *newx;
    SEXP newx3;

    n1 = LENGTH(id2);  /* baseline data set */
    n2 = LENGTH(nid2); /* new data set */

    id	= INTEGER(id2);
    nid = INTEGER(nid2);
    time2 = REAL(time2x);
    ntime = REAL(ntime2);
    x     = REAL(x2);
    indx  = INTEGER(indx2);

    PROTECT(newx3 = duplicate(newx2));
    newx = REAL(newx3);

    /*
    **	There are two indices i= baseline data set, k= additions
    **  For a given subject we might have time intervals of
    **   say (2, 5], [5,9], (12,15], (15,18]
    **   with a newtime, x pairs of (1,105), (15, 202)
    **  In this case the first 3 intervals get a value of 105 and
    **   the last a value of 202.  
    **  The indx variable says where to start for each new addition,
    **   one continues as long as time2 > newtime and id=nid.
    */
    for (k=0; k<n2; k++) {
	for (i=indx[k]-1; i<n1; i++) {
	    if (id[i] != nid[k] || time2[i] <= ntime[k]) break;
	    newx[i] = x[k];
	    }
	}
    
    UNPROTECT(1);
    return(newx3);
    }

/*
** version 2 of the code, used for tdc
**  for each row of the master data (id, time2), return the row of
**  the new data (nid, ntime2) that will provide the new data.
** Based on a last-value-carried forward rule, if the master for an id
**  had time values of 5, 10, 15, 20 and the new data has time values of 7, 15,
**  and 30, the return index would be 0, 1, 1, and 2.  (Covariates change
**  after a time point change.)
*/
	    
SEXP tmerge2(SEXP id2,  SEXP time2x, SEXP nid2, SEXP ntime2) {
    int i, k;
    int n1, n2;
    int oldid;

    int *id, *nid;
    double *time2, 
	   *ntime;
    SEXP index2;
    int  *index;

    n1 = LENGTH(id2);  /* baseline data set */
    n2 = LENGTH(nid2); /* new data set */

    id	= INTEGER(id2);
    nid = INTEGER(nid2);
    time2 = REAL(time2x);
    ntime = REAL(ntime2);

    PROTECT(index2 = allocVector(INTSXP, n1));
    index = INTEGER(index2);

    /*
    ** Every subject in the new data (nid, ntime) will be found in the baseline
    ** data set (id, time2), but not necessarily vice-versa.  
    **   The code walks through the data using i=current row of baseline id,
    ** k = current row of nid.  Visualize each value k as a file folder that 
    ** we want to insert into a file drawer containing the baseline values. 
    ** Folders sorted by time within id, new folder after the baseline if there
    ** is a tied time.
    ** The return vector 'index' is of length n1 (id) and contains the index of
    ** the file folder k that preceeds each i, or -1 if there is no preceeding
    ** insertion with the same id.  We actually return 0 and k+1 for easier
    ** R indexing in the parent.
    **
    **  We cycle between two actions.
    **   0. Increment k 
    **   1. Just found a new value of nid[k].  
    **      a. while (id[i] == prior nid), set index[i]=k and increment i.
    **      b. find the first obs with id[i]==nid[k] and ntime[k] < time[i].
    **     The new obs inserts just in front of 'i' in the file drawer. 
    **     Subjects with this id that preceed the time get a 0.
    **     There might not be such an observation, if the new value would be
    **     the last folder for subject nid[k]; in which case the insertion will
    **     not be used, otherwise set index[i] = k+1.
    **   2. nid[k] = prior value = nid[k-1]
    **      While id[i]==nid[k] and time[i] <= ntime[k], set index[i]=k and
    **      increment i
    **     
    */
    i=0;  /* index for "id" */
    oldid = -1;   /* not anybody */
    for (k=0; k<n2; k++) {
	if (oldid != nid[k]) {
	    for (; i<n1 && (id[i] == oldid); i++) index[i] =k;
	    oldid = nid[k];
	    for(; i<n1 && (id[i] < oldid || 
			   (id[i]== oldid && time2[i] <= ntime[k])); i++)
		index[i] = 0;
	    if (i<n1 && id[i] == oldid) {index[i] = k+1; i++;}
	}	
	else {
	    for (; i<n1 && (id[i] == oldid && time2[i] <= ntime[k]); i++){
		index[i] = k;
	    }
	}
    }	
    for (; i<n1; i++) {
	if (id[i]==oldid) index[i] = k;
	else index[i] =0;
    }

    UNPROTECT(1);
    return(index2);
    }

	    
       
