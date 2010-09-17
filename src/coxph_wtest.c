/* $Id: coxph_wtest.c 11166 2008-11-24 22:10:34Z therneau $ 
** C wrapper for the Cholesky-based Wald test routine
*/
#include "survS.h"
#include "survproto.h"

void coxph_wtest(Sint *nvar2, Sint *ntest, double *var, double *b,
		 double *solve, double *tolerch) {
    int i,j;
    int nvar, df;
    double sum;
    double **var2;
    double *b2;

    nvar = *nvar2;  	/* create a non-pointer version of nvar */
    var2 = dmatrix(var, nvar, nvar); /*make ragged array version of matrix */
    b2 = b;

    cholesky2(var2, nvar, *tolerch);
    df=0;
    for (i=0; i<nvar; i++) if (var2[i][i] >0) df++;  /* count up the df */

    for (i=0; i< *ntest; i++) {
	for (j=0; j<nvar; j++) solve[j] = b[j];
	chsolve2(var2, nvar, solve);   /*solve now has b* var-inverse */

	sum =0;
	for (j=0; j<nvar; j++) sum += b[j]*solve[j];
	b2[i] = sum;                     /* save the result */
	b += nvar;    /*move to next column of b */
	solve += nvar;
        }
    *nvar2 = df;
    }
    
