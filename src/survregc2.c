/*
** Survreg Compute 2 (survregc2)
** Compute the loglik and other parameters for a survreg fit.
**   This is the version for the "callback".  The code below
**   was previously the dolik() routine and it's children found in 
**   survreg3, modified for a joint R/Splus syntax.
** This has exactly the same arguments as survregc2 and returns the
**   same values: see that routine for more comments.  This actually
**   uses expr, rho, and z, however.
*/ 
#include <math.h>
#include "survS.h"
#define SMALL -200  /* exp(-200) is a really small loglik */

double survregc2(int n,          int nvar,     int nstrat,      int whichcase,
		 double *beta,   int dist,     Sint *strat,     double *offset,
		 double *time1,  double *time2, double *status, double *wt,
		 double **covar, double **imat, double **JJ,    double *u, 
		 SEXP expr,      SEXP rho,      double *z,      int nf,
		 Sint *frail,    double *fdiag, double *jdiag ) {

    int person, i,j,k;
    int nvar2;
    int strata;
    double  eta,
	    sigma;
    int     icount;  /* running count of # of interval censored */
    int     fgrp =0; /* the =0 to quiet a compiler warning */
    double  loglik,
	    temp;
    double  temp1, temp2;
    double  sz, zz, zu;
    double  sig2;
    /* add "=0" to keep the compiler from worrying about uninitialized vars */
    /* double g, dg, ddg, dsig, ddsig, dsg; */
    double g=0, dg=0, ddg=0, dsig=0, ddsig=0, dsg=0;
    SEXP   rmat;
    double *funs[5];
    double w;

    nvar2 = nvar + nstrat;
    loglik=0;
    if (whichcase==0) {
	for (i=0; i<nf; i++) {
	    fdiag[i] =0;
	    jdiag[i] =0;
	    }
	for (i=0; i<nvar2+nf; i++) {
	    u[i] =0;
	    for (j=0; j<nvar2; j++) {
		imat[j][i] =0 ;
		JJ[j][i] =0;
		}
	    }
	}

    strata =0;
    sigma = exp(beta[nvar+nf]); /*single scale value case (fixed or estimate)*/
    sig2  = 1/(sigma*sigma);
    loglik =0;

    /*
    ** First, get the array of distribution values
    ** We get them all at once to minimize the S-callback overhead
    */
    icount =n;
    for (person=0; person<n; person++) {
	if (nstrat>1) {
	    strata= strat[person] -1; /*S likes to start counting at 1 */
	    sigma = exp(beta[strata+nvar+nf]);
	    }
	eta =0;
	for (i=0; i<nvar; i++) eta += beta[i] * covar[i][person];
	eta += offset[person];
	if (nf >0){
	    fgrp = frail[person] -1;
	    eta += beta[fgrp]; 
	    }
	z[person] = (time1[person] - eta)/sigma;  
	if (status[person]==3) {
	    z[icount] = (time2[person] - eta)/sigma;
	    icount++;
	    }
	}

    /* 
    ** The result of the eval will be a matrix of 5 rows and n colums, which
    **   we re-index for convenience.  Note that the parent routine has given
    **   us the address of z WITHIN the evaluation frame rho, we just keep
    **   replacing the values it contains; expr then acts like a function of
    **   z.
    ** Actually, if there were any interval censored obs they take up 2 cols;
    **   icount from above contains the actual number of columns used.
    */
    PROTECT(rmat = eval(expr, rho)); 
    funs[0] = REAL(rmat);
    for (i=0; i<4; i++) funs[i+1] = funs[i] + icount;

    /*
    ** calculate the first and second derivative wrt eta,
    **   then the derivatives of the loglik (u, imat, JJ)
    */
    icount =n;
    for (person=0; person<n; person++) {
	if (nstrat>1) {
	    strata= strat[person] -1; /*S likes to start counting at 1 */
	    sigma = exp(beta[strata+nvar]);
	    sig2  = 1/(sigma*sigma);
	    }

	zz = z[person];
	sz = zz * sigma;
	j = status[person];       /*convert to integer */
	switch(j) {
	    case 1:                             /* exact */
		if (funs[2][person] <=0) {
		    /* off the probability scale -- avoid log(0), and set the
		    **  derivatives to gaussian limits (almost any deriv will
		    **  do, since the function value triggers step-halving).
		    */
		    g = SMALL;
		    dg = -zz/sigma;
		    ddg = -1/sigma;
		    dsig =0; ddsig=0; dsg=0;
		    }
		else {
		    g = log(funs[2][person])  - log(sigma);
		    temp1 = funs[3][person]/sigma;
		    temp2 = funs[4][person]*sig2;
		    dg = -temp1;
		    dsig= -(sz*temp1 +1);
		    ddg= temp2 - dg*dg;
		    dsg = sz * temp2 - dg*(1- sz*temp1);
		    ddsig = sz*sz*temp2 + sz*temp1*(1- sz*temp1);
		    }
		break;
	    case 0:                             /* right censored */
		if (funs[1][person] <=0) {
		    g = SMALL;
		    dg = zz/sigma;
		    ddg =0;
		    dsig =0; ddsig=0; dsg=0;
		    }
		else {
		    g = log(funs[1][person]);
		    temp1 = -funs[2][person]/(funs[1][person]*sigma);
		    temp2 = -funs[3][person]*funs[2][person]*sig2/
			               funs[1][person];
		    dg = -temp1;
		    dsig= -sz * temp1;
		    ddg= temp2 - dg*dg;
		    dsg = sz * temp2 - dg*(1+dsig);
		    ddsig = sz*sz*temp2 - dsig*(1+dsig);
		    }
		break;
	    case 2:                             /* left censored */
		if (funs[2][person] <=0) {
		    /* off the probability scale -- avoid log(0) */
		    g = SMALL;
		    dg = -zz/sigma;
		    dsig =0; ddsig=0; dsg=0;
		    ddg =0;
		    }
		else {
		    g = log(funs[0][person]);
		    temp1 = funs[2][person]/(funs[0][person]*sigma);
		    temp2 = funs[3][person]*funs[2][person]*sig2/
			                      funs[0][person];
		    dg= -temp1;
		    dsig= -sz * temp1;
		    ddg= temp2 - dg*dg;
		    dsg = sz * temp2 - dg*(1+dsig);
		    ddsig = sz*sz*temp2 - dsig*(1+dsig);
		    }
		break;
	    case 3:                             /* interval censored */
                zu = z[icount];
                /*stop roundoff in tails*/
		if (zz>0)  temp = funs[1][person] - funs[1][icount]; 
		else       temp = funs[0][icount] - funs[0][person];
		if (temp <=0) {
		    /* off the probability scale -- avoid log(0) */
		    g = SMALL;
		    dg = 1; 
		    ddg =0;
		    dsig =0; ddsig=0; dsg=0;
		    }
		else {
		    funs[3][icount] *= funs[2][icount];  /*f', not f'/f */
		    funs[3][person] *= funs[2][person];
		    g = log(temp);
		    dg  = -(funs[2][icount] -funs[2][person])/(temp*sigma);
		    ddg = (funs[3][icount] -funs[3][person])*sig2/temp - dg*dg;
		    dsig = (zz*funs[2][person] - zu*funs[2][icount])/temp;
		    ddsig= (zu*zu*funs[3][icount] - zz*zz*funs[3][person])
			          /temp - dsig*(1+dsig);
		    dsg = (zu*funs[3][icount] - zz*funs[3][person])/
			       (temp*sigma)  - dg *(1+dsig);
		    }
		icount++;
		break;
	    }
	loglik += g * wt[person];

	/*
	** Now the derivs wrt loglik
	*/
	if (whichcase==1) continue;     /*only needed the loglik */
	w = wt[person];
	if (nf>0) {
	    fgrp = frail[person] -1;
	    u[fgrp] += dg * w;
	    fdiag[fgrp] -= ddg * w;
	    jdiag[fgrp] += dg*dg *w;
	    }
	for (i=0; i<nvar; i++) {
	    temp = dg * covar[i][person] *w;
	    u[i+nf] += temp;
	    for (j=0; j<=i; j++) {
		imat[i][j+nf] -= covar[i][person] *covar[j][person] *ddg *w;
		JJ[i][j+nf]   += temp * covar[j][person] * dg;
		}
	    if (nf>0) {
		imat[i][fgrp] -= covar[i][person] * ddg * w;
		JJ  [i][fgrp] += temp * dg;
		}
	    }

	if (nstrat!=0) {   /* need derivative wrt log sigma */
	    k = strata+nvar;
	    u[k+nf] += w* dsig;
	    for (i=0; i<nvar; i++) {
 		imat[k][i+nf] -= dsg * covar[i][person] * w;
		JJ[k][i+nf]   += dsig* covar[i][person] *dg * w;
		}
	    imat[k][k+nf] -=  ddsig * w;
	    JJ[k][k+nf] += dsig*dsig * w;

	    if (nf>0) {
		imat[k][fgrp] -= dsg * w;
 		JJ  [k][fgrp] += dsig *dg *w;
 		}
	    }
	}

    UNPROTECT(1);   /* release the memory pointed to by funs[] */
    return(loglik);
    }


