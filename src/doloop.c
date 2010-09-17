/* $Id: doloop.c 11357 2009-09-04 15:22:46Z therneau $
**
** Program to mimic a set of nested do loops
**
**  Usual calling sequence would be
**      init_doloop(min,max);
**      while (doloop(nloops, index) >=min) {
**              some calculations
**              }
**
**  The result of this is as though the code had been written for "nloops"
**   nested for loops:
**
**      for (index[0]=min; index[0]<max; index[0]++) {
**          for (index[1]=index[0]+1; index[1]<max; index[1]++) {
**              for (index[2]=index[1]+1; index[2]<max; index[2]++) {
**                    .
**                    .
**
**  min:       The lower limit of the loops
**  max:       The upper limit of the loops
**  nloops:    The number of loops
**  index:   On input, the index of the "last" loop iteration, on
**                  return will contain index of the "next" iteration.
**             (On the first call after init_doloop, input values are ignored)
**
**  returned value:  init_doloop -- none
**
**                   doloop  --     (min-1) signals end of the doloops
**                                  otherwise the value of the innermost loop
**
*/
#include "survS.h"
#include "survproto.h"

static int maxval, minval;
static int firsttime, depth;

void init_doloop(int min, int max) {
    firsttime =1;
    minval = min;
    maxval = max;
    depth =1;
    }

int doloop (int nloops, int *index) {
    register int i;

    if (firsttime ==1) {
	for (i=0; i<nloops; i++)  index[i] =minval+i;
	firsttime =0;
	if (maxval >= (minval+i)) return (minval+i-1);
	    else                  return (minval-1);
	}

    nloops--;
    index[nloops]++;    /*increment the lastmost index */

    if (index[nloops] <= (maxval-depth)) return(index[nloops]);
    else if (nloops ==0)                 return(minval - depth);
	 else {
	    depth++;
	    index[nloops] = doloop(nloops, index) +1;
	    depth--;
	    return(index[nloops]);
	    }
    }
