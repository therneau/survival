/* $Id: survfit4.c 11166 2008-11-24 22:10:34Z therneau $ */
/*
** C routine to do a small computation that is hard in Splus
**
** n = number of observations
** d = number of deaths
** x1, x2 = ingredients in the sums
**
** If d=0, then new x1 = new x2 =1  (fill in value)
**    d=1,      new x1 = 1/x1, 
**	        new x2 = (1/x1)^2
**    d=2,      new x1 = (1/2) [ 1/x1 + 1/(x1 - x2/2)]
**              new x2 = (1/2) [  same terms, squared]
**    d=3       new x1 = (1/3) [ 1/x1 + 1/(x1 - x2/3) + 1/(x1 - 2*x2/3)]
**  etc.
*/

#include "survS.h"

void survfit4(Sint *n,	Sint *dd,  double *x1,  double *x2) {
    double temp, temp1, temp2;
    int i,j;
    double d;

    for (i=0; i< *n; i++) {
	d = dd[i];
	if (d==0) {
	    x1[i] =1;
	    x2[i] =1;
	    }
	else if (d==1){
	    temp = 1/x1[i];
	    x1[i] = temp;
	    x2[i] = temp*temp;
	    }
	else {
	    temp1 = 1/x1[i];
	    temp2 = temp1 * temp1;
	    for (j=1; j<d; j++) {
		temp = 1/(x1[i] - x2[i]*j/d);
		temp1 += temp;
		temp2 += temp*temp;
		}
	    x1[i] = temp1/d;
	    x2[i] = temp2/d;
	    }
	}
    }
