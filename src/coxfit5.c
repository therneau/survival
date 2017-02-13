/* A reentrant version of the Coxfit program, for random effects modeling
**   with reasonable efficiency (I hope).  The important arrays are saved
**   from call to call so as to speed up the process.  The x-matrix itself
**   is the most important of these.
**
** coxfit5_a: Entry and intial iteration step for beta=initial, theta=0
**              (no frailty)
**            Most of the same arguments as coxfit2.
**            Allocate and save arrays in static locations.
** coxfit5_b: Iterate to convergence given an initial value.
** coxfit5_c: Compute residuals and release the saved memory.
**
**     McGilchrist's method for frailty with a fixed theta, but for
**     space savings I assume that many elements of imat are zero
**
**  the input parameters are
**
**       maxiter      :number of iterations
**       nused        :number of people
**       nvar         :number of covariates
**       y[2,n]       :row 1: time of event or censoring for person i
**                    :row 2: status for the ith person    1=dead , 0=censored
**       covar(nv,n)  :covariates for person i.
**                        Note that S sends this in column major order.
**       strata(nstrat):sizes of the strata, cumulative
**       sort 	     :  sort order for the obs, last to first within strata
**       offset(n)    :offset for the linear predictor
**       weights(n)   :case weights
**       eps          :tolerance for convergence.  Iteration continues until
**                       the percent change in loglikelihood is <= eps.
**       tolerch      :tolerance for the Cholesky routines
**       method       : Method 0=Breslow, 1=Efron
**       ptype        : 1 or 3 -- there is a sparse term
**                    : 2 or 3 -- there is a non-sparse term in the model
**       nfrail       : number of frailty groups (sparse terms), 0 if there are
**                        none
**       frail        : a vector containing the frailty groups
**       fbeta        : initial frailty estimates
**       pdiag        : if 0, then for the non-sparse terms only the diagonal
**                        of the variance matrix is penalized, otherwise the
**                        full matrix is used.
**
**  returned parameters
**       means(nv)    : vector of column means of X
**       beta(nv)     : the vector of answers (at start contains initial est)
**       u(nv)        : score vector
**       imat(nv,nv)  : the variance matrix at beta=final
**                      if flag<0, imat is undefined upon return
**       loglik       :loglik at beta=final
**       flag         :success flag  1000  did not converge
**                                   1 to nvar: rank of the solution
**       maxiter      :actual number of iterations used
**       fbeta(nfrail): fitted frailty values
**       fdiag(nfrail + nvar): diagonal of cholesky of the full inverse
**       jmat         : inverse of the cholesky
**       imat         : cholesky of the information matrix
**       expect       : contains the "expected" for each subject
** 
**  work arrays
**       mark(n)
**       wtave(n)
**       score(n)     
**       a(nvar+ nfrail), a2(nvar+nfrail)
**       cmat(nvar,nvar+nfrail)       ragged array
**       cmat2(nvar,nvar+nfrail)
**       fdiag                         the diagonal of the sparse information
**       oldbeta(nvar + nfrail)         always contains the "last iteration"
**
**  the work arrays are passed as a single
**    vector of storage, and then broken out.
**
**  calls functions:  cholesky3, chsolve3, chinv2
**
**  the data must be sorted by ascending time within strata
*/
#include <math.h>
#include <stdio.h>
#include "survS.h"
#include "survproto.h"

static double **covar, **cmat, **cmat2;
static double *mark, *wtave;
static double *a, *oldbeta, *a2;
static double *offset, *weights;
static int    *status, *frail, *sort;
static double *score, *ttime;  /* Hp-UX really doesn't like "time" as a var */
static double *tmean;
static int    ptype, pdiag;
static double *ipen, *upen, logpen;
static Sint   *zflag;

static double **cmatrix(double *, int, int);

void coxfit5_a(Sint *nusedx,     Sint *nvarx,     double *yy, 
	       double *covar2,   double *offset2, double *weights2,
	       Sint *strata,     Sint *sorted,    double *means, 
               double *beta,     double *u,       double *loglik, 
	       Sint *methodx,    Sint *ptype2,    Sint *pdiag2,
	       Sint *nfrail,     Sint *frail2,
               void *fexpr1,     void *fexpr2,    void *rho) {

    int i,j,k, p, istrat;
    int ii; 
    int     nused, nvar;
    int    nf, nvar2;
    
    double  denom, zbeta, risk;
    double  temp, temp2;
    double  ndead;
    double  d2, efron_wt;
    double  method;

    nused = *nusedx;
    nvar  = *nvarx;
    nf= *nfrail;
    method= *methodx;
    nvar2 = nvar + nf;
    ptype = *ptype2;
    pdiag = *pdiag2;

   /*
    **  Allocate storage for the arrays and vectors
    **  Since they will be used later, sizes are based on what will be
    **    needed with the frailty terms.
    */
    if (nvar >0) {
	covar= cmatrix(covar2, nused, nvar);
	cmat = cmatrix(0, nvar2, nvar+1);
	cmat2= cmatrix(0, nvar2, nvar+1);
        }

    a = Calloc(4*nvar2 + 6*nused, double);
    oldbeta = a + nvar2;
    a2 =  oldbeta + nvar2;
    mark = a2 + nvar2;
    wtave= mark + nused;
    weights = wtave+ nused;
    offset  = weights + nused;
    score   = offset + nused;
    tmean   = score + nused;
    ttime    = tmean + nvar2;
    status  = Calloc(2*nused, int);
    sort    = status + nused;
    for (i=0; i<nused; i++) {
	weights[i] = weights2[i];
	offset[i]  = offset2[i];
	status[i]  = yy[nused +i];
	sort[i]    = sorted[i];
	ttime[i]    = yy[i];
        }

    /* scratch space for penalty 
    **    upen needs to be max(nvar, nfrail), 
    **    ipen max(nfrail, nvar(if pdiag=0) or nvar^2 )
    */
    if (nf > nvar) i=nf; else i=nvar;
    if (nf > nvar*nvar) j=nf; else j=nvar*nvar;
    if (pdiag==0)  upen = Calloc(2*i, double);
    else           upen = Calloc(i+j, double);
    ipen = upen + i;
    if (ptype>1)  zflag = Calloc(nvar, Sint);
    else          zflag = Calloc(2, Sint);

    if (nf>0) {
	frail = Calloc(nused, int);
	for (i=0; i<nused; i++) frail[i] = frail2[i];
        }

    /*
    **   Mark(i) contains the number of tied deaths at this point,
    **    for the last person of several tied times. It is zero for
    **    all other points.
    **   Wtave contains the average weight for the deaths
    */
    temp=0;
    j=0;
    istrat=0;
    for (i=0; i<nused; i++)
	mark[i] =0;
    for (i=0; i<nused; ) {
	p = sort[i];
	if (status[p]==1) {
	    temp = 0;
	    ndead=0;
	    for (j=i; j<nused; j++) {
		k = sort[j];
		if ((ttime[k] != ttime[p]) || (j==strata[istrat])) break;
		ndead += status[p];
		temp += weights[k];
		}
	    k=sort[j-1];
	    mark[k] = ndead;
	    wtave[k] = temp/ndead;
	    i=j;
	    }
	else i++;
	if (i==strata[istrat]) istrat++;
	}

    /*
    ** Subtract the mean from each covar, as this makes the regression
    **  much more stable
    */
    for (i=0; i<nvar; i++) {
	temp=0;
	for (p=0; p<nused; p++) temp += covar[i][p];
	temp /= nused;
	means[i] = temp;
	for (p=0; p<nused; p++) covar[i][p] -=temp;
	}

    /*
    ** do the initial iteration step of a no-frailty model
    **   (actually, just a no-sparse-terms model) -- loglik only
    */
    *loglik = 0;
    for (i=0; i<nvar; i++) {
	u[i] =0;
	a[i] = 0;
	a2[i]=0 ;
	}
    denom = 0;
    efron_wt =0;
    istrat=0;
    for (ii=0; ii<nused; ii++) {
	if (ii==strata[istrat]) {
	    denom = 0;
	    for (i=0; i<nvar; i++) a[i] = 0;
	    istrat++;
	    }
	
	p = sort[ii];
	zbeta = offset[p];    /* form the term beta*z   (vector mult) */
	for (i=0; i<nvar; i++)
	    zbeta += beta[i]*covar[i][p];
	zbeta = coxsafe(zbeta);
	risk = exp(zbeta) * weights[p];
	denom += risk;

	for (i=0; i<nvar; i++) a[i] += risk*covar[i][p];
	if (status[p]==1) {
	    efron_wt += risk;
	    *loglik += weights[p]*zbeta;
	    for (i=0; i<nvar; i++) {
	        u[i] += weights[p]*covar[i][p];
		a2[i] +=  risk*covar[i][p];
		}
	    }
	if (mark[p] >0) {  /* once per unique death time */
	    /*
	    ** Trick: when 'method==0' then temp=0, giving Breslow's method
	    */
	    ndead = mark[p];
	    for (k=0; k<ndead; k++) {
		temp = (double)k * method / ndead;
		d2= denom - temp*efron_wt;
		*loglik -= wtave[p] * log(d2);
		for (i=0; i<nvar; i++) {
		    temp2 = (a[i] - temp*a2[i])/ d2;
		    u[i] -= wtave[p] *temp2;
		    }
		}
	    efron_wt =0;
	    for (i=0; i<nvar; i++) a2[i]=0;
	    }
	}   /* end  of accumulation loop */


    /*
    ** add in the penalty terms
    */
    if (ptype==2 || ptype==3) {
	/* there are non-sparse terms */
	cox_callback(2, beta, upen, ipen, &logpen, zflag, nvar, fexpr2, rho);
	*loglik += logpen;
        }
    }

/********************************************************************/

/*
** This call is used for iteration
*/
void coxfit5_b(Sint *maxiter, Sint *nusedx, Sint *nvarx, 
	       Sint *strata, double *beta, double *u,
	       double *imat2,  double *jmat2, double *loglik, 
	       Sint *flag,  double *eps, double *tolerch, Sint *methodx, 
	       Sint *nfrail, double *fbeta, double *fdiag,
               void *fexpr1, void *fexpr2, void *rho)
{

    int i,j,k, p;
    int ii, istrat, ip;
    int     iter;
    int     nused, nvar;
    int    nf, nvar2;
    int    fgrp=0;
    int    halving;

    double  denom=0, zbeta, risk;
    double  temp, temp2;
    double  newlk=0;
    double  d2, efron_wt=0;
    double  method;
    double  **imat, **jmat;
    double  ndead;

    nused = *nusedx;
    nvar  = *nvarx;
    nf= *nfrail;
    method= *methodx;
    nvar2 = nvar + nf;
    if (nvar >0) {
        imat  = dmatrix(imat2, nvar2, nvar);
	jmat  = dmatrix(jmat2, nvar2, nvar);
        }
    else {
	imat = 0;   /*never used, but passed as dummy to chol */
	jmat = 0;
        }

    for (i=0; i<nf; i++) oldbeta[i] = fbeta[i];
    for (i=0; i<nvar; i++) oldbeta[i+nf] = beta[i];

    halving =0 ;             /* =1 when in the midst of "step halving" */
    for (iter=0; iter<=*maxiter; iter++) {
	newlk = 0;
	for (i=0; i<nf; i++) fdiag[i] =0;
	for (i=0; i<nvar2; i++) {
	    u[i] =0;
	    for (j=0; j<nvar; j++)
		    jmat[j][i] =0 ;
            }

	istrat=0;
	for (ip=0; ip<nused; ip++) {
	    p = sort[ip];
	    if (ip==0 || ip==strata[istrat]) {
		efron_wt =0;
		denom = 0;
		for (i=0; i<nvar2; i++) {
		    a[i] = 0;
		    a2[i]=0 ;
		    for (j=0; j<nvar; j++) {
			cmat[j][i] = 0;
			cmat2[j][i]= 0;
                        }
		    }
		}
	    if (ip==strata[istrat]) istrat++;

	    if (nf>0) {
		fgrp = frail[p] -1;
		zbeta = offset[p] + fbeta[fgrp];
	        }
	    else zbeta = offset[p];

	    for (i=0; i<nvar; i++)
		zbeta += beta[i]*covar[i][p];
	    zbeta = coxsafe(zbeta);
	    score[p] = exp(zbeta);
	    risk = score[p] * weights[p];
	    denom += risk;

	    if (nf>0) a[fgrp] += risk;
	    for (i=0; i<nvar; i++) {
		a[i+nf] += risk*covar[i][p];
		if (nf>0) cmat[i][fgrp] += risk*covar[i][p];
		for (j=0; j<=i; j++)
		    cmat[i][j+nf] += risk*covar[i][p]*covar[j][p];
		}

	    if (status[p]==1) {
		efron_wt += risk;
		newlk += weights[p] *zbeta;
		if (nf>0) {
		    u[fgrp] += weights[p];
		    a2[fgrp] += risk;
		    }
		for (i=0; i<nvar; i++) {
		    u[i+nf] += weights[p] *covar[i][p];
		    a2[i+nf] +=  risk*covar[i][p];
		    if (nf>0) cmat2[i][fgrp] += risk*covar[i][p];	
		    for (j=0; j<=i; j++)
			cmat2[i][j+nf] += risk*covar[i][p]*covar[j][p];
   		    }
		}

	    if (mark[p] >0) {  /* once per unique death time */
		ndead = mark[p];
		for (k=0; k<ndead; k++) {
		    temp = (double)k* method / ndead;
		    d2= denom - temp*efron_wt;
		    newlk -= wtave[p] *log(d2);

		    for (i=0; i<nvar2; i++) {  /* by row of full matrix */
			temp2 = (a[i] - temp*a2[i])/d2;
			tmean[i] = temp2;
			u[i] -= wtave[p] *temp2;
			if (i<nf) fdiag[i] += temp2 * (1-temp2);
			else {
			    ii = i-nf;     /*actual row in c/j storage space */
			    for (j=0; j<=i; j++) 
				jmat[ii][j] +=  wtave[p]*(
                                   (cmat[ii][j] - temp*cmat2[ii][j]) /d2 -
                                          temp2*tmean[j]);
			    }
                        }
		    }
		efron_wt =0;
		for (i=0; i<nvar2; i++) {
		    a2[i]=0;
		    for (j=0; j<nvar; j++)  cmat2[j][i]=0;
		    }
		}
	    }   /* end  of accumulation loop  */

	/*
	** Add in the penalty terms
        */
	if (ptype==1 || ptype==3) {
	    /* there are sparse terms */
    	    cox_callback(1, fbeta, upen, ipen, &logpen, zflag, nf, fexpr1,rho); 
	    if (zflag[0] ==1) {  /* force terms to zero */
		for (i=0; i<nf; i++) {
		    u[i]=0;
		    fdiag[i] =1;
		    for (j=0; j<nvar; j++) jmat[j][i]=0;
		    }
	        }
	    else {
		for (i=0; i<nf; i++) {
		    u[i] += upen[i];
		    fdiag[i] += ipen[i];
		    }
		newlk += logpen;
	        }
	    }

	if (ptype==2 || ptype==3) {
	    /* there are non-sparse terms */
	    cox_callback(2, beta, upen, ipen, &logpen, zflag, nvar, fexpr2,rho);
	    newlk += logpen;
	    if (pdiag==0) {
		for (i=0; i<nvar; i++) {
		    u[i+nf] += upen[i];
		    jmat[i][i+nf] += ipen[i];
		    }
	        }
	    else {
		k =0;
		for (i=0; i<nvar; i++) {
		    u[i+nf] += upen[i];
		    for (j=nf; j<nvar2; j++) jmat[i][j] += ipen[k++];
		    }
	        }
	    for (i=0; i<nvar; i++) {
		if (zflag[i] ==1) {
		    u[i+nf]=0;
		    for (j=0; j<i; j++) jmat[i][j+nf]=0;
		    jmat[i+nf][i] =1;
		    }
	        }
	    }

	/* am I done?
	**   update the betas and test for convergence
	*/
	*flag = cholesky3(jmat, nvar2, nf, fdiag, *tolerch);
	if (fabs(1-(*loglik/newlk))<=*eps && halving==0) { /* all done */
	    *loglik = newlk;
	    for (i=0; i<nvar; i++) {
	        for (j=0; j<nvar2; j++)  imat[i][j] = jmat[i][j];
	        }
	    chinv3(jmat, nvar2, nf, fdiag);
	    for (i=nf; i<nvar2; i++) {       /*nicer output for S user */
	        fdiag[i] = jmat[i-nf][i];
	        jmat[i-nf][i] =1;
		imat[i-nf][i] =1;
		for (j=i+1; j<nvar2; j++) {
		    jmat[i-nf][j] = 0;
		    imat[i-nf][j] = 0;
		    }
	        }

	    *maxiter = iter;
	    return;
	    }

	if (iter==*maxiter) break;  /*skip the step halving and etc */

	if (iter>0 && newlk < *loglik)   {    /*it is not converging ! */
		halving =1;
		for (i=0; i<nvar; i++)
		    beta[i] = (oldbeta[i+nf] + beta[i]) /2; 
		for (i=0; i<nf; i++)
		    fbeta[i] = (oldbeta[i] + fbeta[i])/2;
		}
	    else {
		halving=0;
		*loglik = newlk;
		chsolve3(jmat,nvar2, nf, fdiag, u);

		j=0;
		for (i=0; i<nvar; i++) {
		    oldbeta[i+nf] = beta[i];
		    beta[i] += u[i+nf];
		    }
		for (i=0; i<nf; i++) {
		    oldbeta[i] = fbeta[i];
		    fbeta[i] += u[i];
		    }
		}
	}   /* return for another iteration */

    *loglik = newlk;
    for (i=0; i<nvar; i++) 
	for (j=0; j<nvar2; j++) {
            imat[i][j] = jmat[i][j];
        }
    chinv3(jmat, nvar2, nf, fdiag);
    for (i=nf; i<nvar2; i++) {       /*nicer output for S user */
	fdiag[i] = jmat[i-nf][i];
	jmat[i-nf][i] =1;
  	imat[i-nf][i] =1;
  	for (j=i+1; j<nvar2; j++) {
	    jmat[i-nf][j] = 0;
  	    imat[i-nf][j] = 0;
  	    }
        }
    *flag= 1000;
    return;
    }


static double **cmatrix(double *data, int ncol, int nrow)
    {

    int i,j;
    double **pointer;
    double *temp;
 
    pointer = Calloc(nrow, double *);
    temp = Calloc(nrow*ncol, double);
    if (data==0){
	for (i=0; i<nrow; i++) {
	    pointer[i] = temp;
	    temp += ncol;
            }
        }
    else {
	for (i=0; i<nrow; i++) {
	    pointer[i] = temp;
	    for (j=0; j<ncol; j++) *temp++ = *data++;
	    }
        }
    return(pointer);
    }

static void cmatrix_free(double **data) 
{
    Free(*data);
    Free(data);
    }


void coxfit5_c (Sint *nusedx, Sint *nvar, Sint *strata, Sint *methodx, 
		double *expect) {
    double hazard, 
           denom,
           temp, temp2,
           efron_wt,
           ndead, 
           hazard2;
    int    p,
           nused,
           method,
	   ip, istrat,
           i, j;

    nused = *nusedx;
    method= *methodx;

    /*
    ** compute the expected number of events for each subject 
    */
    istrat=0;
    denom =0;
    for (ip=0; ip<nused; ip++) {
	p = sort[ip];

        if (ip==strata[istrat]) {
	    denom=0;
	    istrat++;
	    }
        denom += score[p] * weights[p];

	if (mark[p] >0) {
	    /*
	    ** Compute the size of the hazard jump at this point, with the
	    **  total jump saved (temporarily) in "expect", and the Efron
	    **  amount in "weights".  It applies to deaths at this point.
	    */
	    ndead = mark[p];
	    temp2  = 0; 
	    efron_wt =0;
	    for (j=0; j<ndead; j++) {  
		/* walk backwards in ip, over the tied deaths */
		i = sort[ip-j];
		efron_wt += score[i]* weights[i];
		temp2    += weights[i];
	        }
        
            if (ndead<2 || method==0)  {
		expect[p]  = temp2/denom;
		weights[p] = temp2/denom;
		}
	    else {
		hazard =0;
		hazard2=0;
		temp2 /= ndead;
                for (j=0; j<ndead; j++) {
                    temp = j /ndead;
                    hazard +=  temp2/(denom - efron_wt* temp);
                    hazard2+=  temp2*(1-temp)/(denom - efron_wt* temp);
                    }
		expect[p]  = hazard;
		weights[p] = hazard2;
	        }
	    }
	}
    /*
    ** Now compute cumulative hazard,
    **  and store in the "expect" vector
    */
    hazard=0;
    for (ip=nused-1; ip>=0; ) {
	p = sort[ip];
	if (status[p] >0) {
	    ndead = mark[p];
	    temp  =  expect[p];  
	    hazard2 =weights[p];

	    for (j=0; j<ndead; j++) {
		i = sort[ip-j];
		expect[i] = score[i]*(hazard + hazard2);
		}
	    ip -= ndead;
	    hazard += temp;
	    }
	else {
	    expect[p] = hazard * score[p];
	    ip--;
	    }

	if (strata[istrat]==ip) {
	    hazard=0;
	    istrat--;
	    }
	}
    
    /*
    ** Free up the extra memory
    */
    Free(zflag);
    Free(upen);
    Free(status);
    Free(a);
    if (*nvar > 0) {
	cmatrix_free(cmat2);
	cmatrix_free(cmat);
	cmatrix_free(covar);
        }
    }





