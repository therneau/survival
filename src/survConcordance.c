/*
**    $Id: survConcordance.c 11166 2008-11-24 22:10:34Z therneau $
**
**  For each observation, we want to know, for the subset of observations
**     with longer survival (and only those)
**          number with smaller, bigger, and tied x values
**
**  The input data is sorted, largest survival to smallest survival
**
**     n	number of time/status/x values
**     time
**     status   needed to keep track of tied survival times
**     x        vector of scores
**     n2       number of unique x values
**     x2       sorted vector of unique x values, smallest to largest
**
**     temp     scratch vector of length 2* n2
**
**  returned
**     result   number concordant, discordant, tied survival, tied x but
**                not tied survival, and incomparable times
**     (bigger survival + smaller risk score = concordant)
*/
#include "survS.h"
#include <stdio.h>
void survConcordance(Sint *np,    double *time,  Sint *status, 
		     double *x,   Sint *n2p,     double *x2,
		     Sint   *temp,Sint *result) {
    int i, j, k=0;
    int start, end;
    int n, n2;
    Sint *count1, *count2, *count;
    int tdeath; 
    int nright, nsame;

    n = *np;
    n2= *n2p;
    count1 = &(temp[0]);
    count2 = &(temp[n2]);
    for (i=0; i<5; i++) result[i] =0;  /* redundant I think */
    for (i=0; i<n2; i++) count1[i]=0;

    /* 
    ** The heart of the algorithm is to think of the ordered list of
    **   unique values as a balanced binary tree.  (Credit to Brad Broom
    **   of Rice U for this idea).
    **  For any node k, below it to the left all values are < x2[k],
    **    and below to the right all values are > x2[k].  (Draw a picture).
    **  The root of the tree is element k= floor((n2-1)/2), with value x2[k].
    **  In general, for any subtree that "owns" elements i to j, the root
    **    of that subtree is element k= floor((i+j)/2), whose left subtree
    **    owns elements i to k-1 of the tree, and right subtree owns elements
    **    k+1 to j.
    **     
    **  As we update, count[i] will be the number of data values in this
    **     node and all nodes below.
    **
    **  We walk through the data one survival time at at time, comparing each
    **    to all the survival times above it.
    **  If the time is censored, all those above are "incomparable".
    **  Otherwise, we need to find the position of x[i], among x[1: (i-1)]
    **  We do this by updating the counts in the binary tree.  The count
    **  vector contains the number of x[0 to i] that are in or below any
    **  given node k of the binary tree.
    ** 
    **  Tied death times are a nuisance; we have to refrain from updating
    **   the counts until the end of each set of them.  Thus a vector
    **   count1 (up to date) and count2 (lagged).
    **  nright = sum(# values to the right, each time I take a left branch)
    */

    tdeath =0;     /* current count of tied deaths */
    for (i=0; i<n; i++) {
	if (status[i] > 0) {
	    /*
	    ** Walk the tree a first time, to count this observation's
	    **   position
	    */
	    nright = 0;
	    start = 0; end= n2-1;  /*start to end of sublist being looked at */
	    if (tdeath==0) count=count1;   /* use the appropriate count */
	    else           count=count2;
	    while(start <= end) {
		k = (start+end)/2;
		if (x[i] == x2[k]) break;
		if (x[i] < x2[k])  {
		    /* take the left branch (smaller numbers) */
		    end = k-1;
		    nright = nright + (count[k] - count[(start+end)/2]);
		    }
		else  start = k+1;  /*right branch */
                }

	    /*
	    ** At this point x[i] = x2[k]; we've found the number in the 
	    **  x2 list
	    */
	    nsame = count[k];   /*provisional */
	    if (k<end) {  /* there is a right hand branch below this node*/
		j = count[(k+1+end)/2];  /*number to the right */
		nsame = nsame -j;
		nright= nright+j;
                }
	    if (k > start)  /* there is a left hand branch below here */
		nsame = nsame - count[(start+k-1)/2];

	    result[3] += nsame;
            result[1] += nright;  /* # values bigger than x[i] */
            result[0] +=  i - (tdeath + nsame + nright); /* # smaller */

            /* Is the next survival time tied with this one? */
            if (i<(n-1) && status[i+1]>0 &&(time[i] == time[i+1])) {
                tdeath += 1;  /* Yes it is */
		if (tdeath==1) {
		    for (j=0; j<n2; j++) count2[j] = count1[j];
		    }
		}
            else {
		result[2] += (tdeath * (tdeath+1))/2;
                tdeath =0;
                }
            }
        else { 
	    /* 
	    ** censored survival time
	    ** All those above it on the list are "incomparable"
	    */
	    tdeath =0;
	    result[4] += i;
            }

	/*
	** Now, walk the list one more time, updating count1
	*/
	start = 0; end= n2-1;  /*start to end of sublist being looked at */
	while(start <= end) {
	    k = (start+end)/2;
	    count1[k]++;
	    if (x[i] == x2[k]) break;
	    if (x[i] < x2[k]) 	end = k-1;  /* left branch */
	    else  start = k+1;  /*right branch */
	    }
        }
    }

