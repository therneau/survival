/*
** HOw the routines fit together
**  survreg.fit, the S function, calls either 
**     survreg6 -- no penalized terms
**     survreg7 -- with penalized terms
**  Each of those routines calls either
**     survregc1 -- one of the 3 built in distributions
**     survregc2 -- a user-written distribution
**   to do the loglike compuations.  
**
**  In addition, survreg7 calls survpenal to deal with the call-back for
**    penalized terms.
**
** Survreg Compute 1 (survregc1)
** Compute the loglik and other parameters for a survreg fit.
**   This is the version for the 3 built-in distributions.  The code below
**   was previously the dolik() routine and it's children found in 
**   survreg2/survreg4.
** The expr, rho and dummy arguments are unused in this routine.
**
**  If there is a sparse term then nf = #levels for the term (there is 
**    never more than 1 sparse term).  The frail variable then contains the
**    group number for each subject, columns 1-nf of the beta are the
**    coefficients for the sparse term (1 per level, no contrasts!).
**
**  The beta vector will have nf + nvar + nstrat coefficients, in that order.
**    The last set are the estimate(s) of variance.  For a fixed variance model
**    nstrat =0, by far the most common case is nvar=1.  When nvar=0, the
**    value for the fixed variance is tacked onto beta, giving it an intial
**    length one longer than one might expect.
**  nvar = # variables in the X matrix (covar).  Variables imat and
**    JJ contain the dense portion of the information and J'J matrices, and the
**    sparse portions of the matrices are in fdiag and jdiag. So imat will
**    have dimension (nvar + nstrat) by (nf + nvar + nstrat).
**
**  When whichcase==1, the parent routine only needs a loglik, which saves some
**    computation.
**
**  See the technical report for all the long formulas that define the
**    derivatives.  (Mayo Biostatistics, TR #53)
*/ 
#include "survS.h"
#include <math.h>
#define SMALL -200 /* what to use for log(f(x)) when f(x) gives a zero,
		       i.e., the calling made a really bad guess for beta */

static void exvalue_d(double z, double ans[4], int j);
static void logistic_d(double z, double ans[4], int j);
static void gauss_d(double z, double ans[4], int j);
static void (*sreg_gg)();

#define  SPI    2.506628274631001     /* sqrt(2*pi) */
#define  ROOT_2 1.414213562373095

double survregc1(int n,          int nvar,      int nstrat,     int whichcase,
		 double *beta,   int dist,      Sint *strat,    double *offset,
		 double *time1,  double *time2, double *status, double *wt,
		 double **covar, double **imat, double **JJ,    double *u, 
		 SEXP expr,      SEXP rho,      double *dummy,  int nf,
		 Sint *frail,    double *fdiag, double *jdiag ) {
    
    int person, i,j,k;
    int nvar2;        /* nvar + nstrat */
    int nvar3;        /* nvar2 + nf */
    int strata;
    double  eta,
	    sigma;
    double  z, zu,
	    loglik,
	    temp, temp2;
    double  sz;
    double  sig2;
    double  funs[4], ufun[4];
    int     fgrp =0;        /* the =0 to quiet a compiler warning */
    double  w;
    /* add "=0" to keep the compiler from worrying about uninitialized vars */
    double g=0, dg=0, ddg=0, dsig=0, ddsig=0, dsg=0;

    switch(dist) {
	case 1: sreg_gg = exvalue_d;  break;
	case 2: sreg_gg = logistic_d; break;
	case 3: sreg_gg = gauss_d;    break;
	}

    nvar2 = nvar + nstrat;
    nvar3 = nvar2 + nf;

    loglik =0;
    if (whichcase==0) {
	for (i=0; i<nf; i++) {
	    fdiag[i] =0;
	    jdiag[i] =0;
	    }
	for (i=0; i<nvar3; i++) {
	    u[i] =0;
	    for (j=0; j<nvar2; j++) {
		imat[j][i] =0 ;
		JJ[j][i] =0;
		}
	    }
	}

    /*
    ** calculate the first and second derivative wrt eta,
    **   then the derivatives of the loglik (u, imat, JJ)
    */
    strata =0;
    sigma = exp(beta[nvar+nf]); 
    sig2  = 1/(sigma*sigma);
    for (person=0; person<n; person++) {
	if (nstrat>1) {
	    /* 
	    ** multiple scales: pick the right sigma for this obs
	    ** The more common case of a single scale is set 6 lines above
	    */
	    strata= strat[person] -1; /*S likes to start counting at 1 */
	    sigma = exp(beta[strata+nvar+nf]);
	    sig2  = 1/(sigma*sigma);
	    }
	eta =0;
	for (i=0; i<nvar; i++) eta += beta[i+nf] * covar[i][person];
	eta += offset[person];
	if (nf >0){
	    fgrp = frail[person] -1;
	    eta += beta[fgrp]; 
	    }
	sz = (time1[person] - eta);  /*   sigma * z  */
	z = sz /sigma;

	j = status[person];       /*convert to integer */
	switch(j) {
	    case 1:                             /* exact */
		(*sreg_gg)(z, funs,1);
		if (funs[1] <=0) {
		    /* off the probability scale -- avoid log(0), and set the
		    **  derivatives to gaussian limits (almost any deriv will
		    **  do, since the function value triggers step-halving).
		    */
		    g = SMALL;
		    dg = -z/sigma;
		    ddg = -1/sigma;
		    dsig =0; ddsig=0; dsg=0;
		    }
		else {
		    g = log(funs[1])  - log(sigma);
		    temp = funs[2]/sigma;
		    temp2= funs[3]*sig2;
		    dg = -temp;
		    dsig= -temp*sz;
		    ddg= temp2 - dg*dg;
		    dsg = sz * temp2 - dg*(dsig +1);
		    ddsig = sz*sz* temp2 - dsig*(1+dsig);
		    dsig -= 1;
		    }
		break;
	    case 0:                             /* right censored */
		(*sreg_gg)(z, funs,2);
		if (funs[1] <=0) {
		    g = SMALL;
		    dg = z/sigma;
		    ddg =0;
		    dsig =0; ddsig=0; dsg=0;
		    }
		else {
		    g = log(funs[1]);
		    temp = -funs[2]/(funs[1]*sigma);
		    temp2= -funs[3]*sig2/funs[1];
		    dg = -temp;
		    dsig= -temp*sz;
		    ddg= temp2 - dg*dg;
		    dsg = sz * temp2 - dg*(dsig +1);
		    ddsig = sz*sz* temp2 - dsig*(1+dsig);
		    }
		break;
	    case 2:                             /* left censored */
		(*sreg_gg)(z, funs,2);
		if (funs[0] <=0) {
		    /* off the probability scale -- avoid log(0) */
		    g = SMALL;
		    dg = -z/sigma;
		    dsig =0; ddsig=0; dsg=0;
		    ddg =0;
		    }
		else {
		    g = log(funs[0]);
		    temp = funs[2]/(funs[0]*sigma);
		    temp2= funs[3]*sig2/funs[0];
		    dg = -temp;
		    dsig= -temp*sz;
		    ddg= temp2 - dg*dg;
		    dsg = sz * temp2 - dg*(dsig +1);
		    ddsig = sz*sz* temp2 - dsig*(1+dsig);
		    }
		break;
	    case 3:                             /* interval censored */
		zu = (time2[person] - eta)/sigma;  /*upper endpoint */
		(*sreg_gg)(z, funs, 2);
		(*sreg_gg)(zu,ufun ,2);
		if (z>0)  temp = funs[1] - ufun[1]; /*stop roundoff in tails*/
		else      temp = ufun[0] - funs[0];
		if (temp <=0) {
		    /* off the probability scale -- avoid log(0) */
		    g = SMALL;
		    dg = 1; 
		    ddg =0;
		    dsig =0; ddsig=0; dsg=0;
		    }
		else {
		    g = log(temp);
		    dg  = -(ufun[2] - funs[2])/(temp*sigma);
		    ddg = (ufun[3] - funs[3])*sig2/temp - dg*dg;
		    dsig = (z*funs[2] - zu*ufun[2])/temp;
		    ddsig= ((zu*zu*ufun[3] - z*z*funs[3])/temp) -
			                dsig*(1+dsig);
		    dsg = ((zu*ufun[3] - z*funs[3])/ (temp*sigma)) -
				      dg * (dsig +1);
		    }
		break;
	    }
	loglik += g * wt[person];
/*if (person<8) fprintf(stderr, "i=%d, g=%g, dg=%g, ddg=%g, dsg=%g\n",
		 person, g, dg, ddg, dsg);*/
	/*
	** Now the derivs wrt loglik
	**  Remember that the "x" for a sparse term is 1
	*/
	if (whichcase==1) continue;     /*only needed the loglik */
	w = wt[person];
	if (nf>0) {
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

    return(loglik);
    }


/*
**  Case      ans[0]    ans[1]       ans[2]     ans[3]
**   1                    f          f'/f        f''/ f
**   2          F        1-F         f           f'
**
**  We do both F and 1-F to avoid the error in (1-F) for F near 1
*/

static void logistic_d(double z, double ans[4], int j)
    {
    double w, temp;
    int    sign, ii;

    /*
    ** The symmetry of the logistic allows me to be careful, and never take
    **  exp(large number).  This routine should be very accurate.
    */
    if (z>0)  {
	w = exp(-z);
	sign = -1;
	ii=0;
	}
    else {
	w = exp(z);
	sign = 1;
	ii=1;
	}
    temp = 1+w;
    switch(j) {
	case 1:  ans[1] = w/(temp*temp);
		 ans[2] = sign*(1-w)/temp;
		 ans[3] = (w*w -4*w +1)/(temp*temp);
		 break;
	case 2:  ans[1-ii] = w/temp;
		 ans[ii]   = 1/temp;
		 ans[2] = w/(temp*temp);
		 ans[3] = sign*ans[2]*(1-w)/temp;
		 break;
	}
    }

static void gauss_d(double z, double ans[4], int j)
    {
    double f;

    f = exp(-z*z/2) /SPI;
    switch(j) {
	case 1: ans[1] =f;
		ans[2] = -z;
		ans[3] = z*z -1;
		break;
	case 2: if (z>0) {
		    ans[0] = (1 + erf(z/ROOT_2))/2;
		    ans[1] =  erfc(z/ROOT_2) /2;
		    }
		else {
		    ans[1] = (1 + erf(-z/ROOT_2))/2;
		    ans[0] =  erfc(-z/ROOT_2) /2;
		    }
		ans[2] = f;
		ans[3] = -z*f;
		break;
	}
    }

/*
** In the Gaussian and logistic cases, I could avoid numeric disaster by only
**   evaluating exp(x) for x<0.  By symmetry, I got what I need for
**   x >0.  The extreme value dist, howerver, is asymmetric, and I don't yet 
**   see the appropriate numeric tricks.
** Perhaps a Taylor series will could be used for large z.
*/

static void exvalue_d(double z, double ans[4], int j)
    {
    double temp;
    double w;
    if (z < SMALL) w= exp(SMALL);
    else if (-z < SMALL) w = exp(-SMALL);  /* stop infinite answers */
    else   w = exp(z);

    temp = exp(-w);
    switch(j) {
	case 1:  ans[1] = w*temp;
		 ans[2] = 1-w;
		 ans[3] = w*(w-3) +1;
		 break;
	case 2:  ans[0] = 1-temp;
		 ans[1] = temp;
		 ans[2] = w*temp;
		 ans[3] = w*temp*(1-w);
		 break;
	}
    }

