/* $Id: agexact.c 11166 2008-11-24 22:10:34Z therneau $ */
/*
** Anderson-Gill formulation of the cox Model
**   Do an exact calculation of the partial likelihood. (CPU city!)
**
**  the input parameters are
**
**       maxiter      :number of iterations
**       nused        :number of people
**       nvar         :number of covariates
**       start(n)     :each row covers the time interval (start,stop]
**       stop(n)      :
**       event(n)     :was there an event at 'stop':1=dead , 0=censored
**       covar(nv,n)  :covariates for person i.
**                        Note that S sends this in column major order.
**       strata(n)    :marks the strata.  Will be 1 if this person is the
**                       last one in a strata.  If there are no strata, the
**                       vector can be identically zero, since the nth person's
**                       value is always assumed to be = to 1.
**       offset(n)    :linear offset
**       eps          :tolerance for convergence.  Iteration continues until
**                       the percent change in loglikelihood is <= eps.
**       tol_chol     : tolerance for the Cholesky routine
**
**  returned parameters
**       means(nv)    :column means of the X matrix
**       beta(nv)     :the vector of answers (at start contains initial est)
**       u            :the first derivative vector at solution
**       imat(nv,nv)  :the variance matrix at beta=final, also a ragged array
**                      if flag<0, imat is undefined upon return
**       loglik(2)    :loglik at beta=initial values, at beta=final
**       sctest       :the score test at beta=initial
**       flag         :success flag  1000  did not converge
**                                   1 to nvar: rank of the solution
**       maxiter      :actual number of iterations used
**
**  work arrays
**       score(n)              the score exp(beta*z)
**       a(nvar)
**       cmat(nvar,nvar)       ragged array
**       newbeta(nvar)         always contains the "next iteration"
**
**  the 4 arrays score, a, cmat, and newbeta are passed as a single
**    vector of storage, and then broken out.
**
**  calls functions:  cholesky2, chsolve2, chinv2
**
**  the data must be sorted by ascending time within strata, deaths before
**          living within tied times.
*/
#include <math.h>
#include "survS.h"
#include "survproto.h"

void agexact(Sint *maxiter,  Sint *nusedx,   Sint *nvarx,   double *start, 
	     double *stop,   Sint *event,    double *covar2,double *offset, 
	     Sint   *strata, double *means,  double *beta,  double *u, 
	     double *imat2,  double loglik[2], Sint *flag,  double *work, 
	     Sint   *work2,  double *eps,    double *tol_chol, double *sctest)
{
    int i,j,k, l, person;
    int     iter;
    int     n, nvar;

    double **covar, **cmat, **imat;  /*ragged array versions*/
    double *a, *newbeta;
    double *score, *newvar;
    double  denom, zbeta, weight;
    double  time;
    double  temp;
    double  newlk =0;
    int     halving;    /*are we doing step halving at the moment? */
    int     nrisk, deaths;
    int *index, *atrisk;

    n = *nusedx;
    nvar  = *nvarx;
    /*
    **  Set up the ragged arrays
    */
    covar= dmatrix(covar2, n, nvar);
    imat = dmatrix(imat2,  nvar, nvar);
    cmat = dmatrix(work,   nvar, nvar);
    a = work + nvar*nvar;
    newbeta = a + nvar;
    score   = newbeta + nvar;
    newvar  = score + n;
    index =  (int *) work2;
    atrisk=  index+n;

    /*
    ** Subtract the mean from each covar, as this makes the regression
    **  much more stable
    */
    for (i=0; i<nvar; i++) {
	temp=0;
	for (person=0; person<n; person++) temp += covar[i][person];
	temp /= n;
	means[i] = temp;
	for (person=0; person<n; person++) covar[i][person] -= temp;
	}

    /*
    ** do the initial iteration step
    */
    loglik[1] =0;
    for (i=0; i<nvar; i++) {
	u[i] =0;
	for (j=0; j<nvar; j++)
	    imat[i][j] =0 ;
	}

    for (person=0; person<n; person++) {
	zbeta = 0;      /* form the term beta*z   (vector mult) */
	for (i=0; i<nvar; i++)
	    zbeta += beta[i]*covar[i][person];
	score[person] = exp(zbeta + offset[person]);
        }

    for (person=0; person<n;) {
	if (event[person]==0) person++;
	else {
	    denom =0;
	    for (i=0; i<nvar; i++) {
		a[i] =0;
		for (j=0; j<nvar; j++)  cmat[i][j]=0;
		}

	    /*
	    ** Compute whom is in the risk group, and #deaths
	    */
	    nrisk=0;
	    deaths=0;
	    time = stop[person];
	    for (k=person; k<n; k++) {
		if (stop[k]==time) deaths += event[k];
		if (start[k] < time) {
		    atrisk[nrisk]=k;
		    nrisk++;
		    }
		if (strata[k]==1) break;
		}

	    /*
	    ** compute the mean and covariance over the risk set (a and c)
	    **   It's fast if #deaths=1
	    */
	    if (deaths==1) {
		for (l=0; l<nrisk; l++) {
		    k = atrisk[l];
		    weight = score[k];
		    denom += weight;
		    for (i=0; i<nvar; i++) {
			a[i] = a[i] + weight*covar[i][k];
			for (j=0; j<=i; j++)
			    cmat[i][j] += weight*covar[i][k]*covar[j][k];
			}
		     }
		}
	    else {
		/*
		** for each unique subset of size "deaths" from the risk set,
		**  the "new variable" is the sum of x over that set.  It's
		**  weight is the product of the weights.  I want a to contain
		**  the weighted sum and c the weighted ss of this new var.
		*/
		init_doloop(0,nrisk);
		while(doloop(deaths, index) >=0) {
		    for (i=0; i<nvar; i++) newvar[i]=0;
		    weight =1;
		    for (l=0; l<deaths; l++) {
			k = atrisk[index[l]];
			weight *= score[k];
			for (i=0; i<nvar; i++)  newvar[i]+= covar[i][k];
			}
		    denom += weight;
		    for (i=0; i<nvar; i++) {
			a[i] = a[i] + weight*newvar[i];
			for (j=0; j<=i; j++)
			    cmat[i][j] += weight*newvar[i]*newvar[j];
			}
		     }
		}
	    /*
	    ** Add results into u and imat
	    */
	    loglik[1] -= log(denom);
	    for (i=0; i<nvar; i++) {
		u[i] -=  a[i]/denom;
		for (j=0; j<=i; j++)
		    imat[j][i] += (cmat[i][j] - a[i]*a[j]/denom)/denom;
		}
	    for (k=person; k<n && stop[k]==time; k++) {
		if (event[k]==1) {
		    loglik[1] +=  log(score[k]);
		    for (i=0; i<nvar; i++) u[i] += covar[i][k];
		    }
		person++;
		if (strata[k]==1) break;
		}
	    }
	}   /* end  of accumulation loop */

    loglik[0] = loglik[1];   /* save the loglik for iteration zero  */

    /* am I done?
    **   update the betas and test for convergence
    */

    for (i=0; i<nvar; i++) /*use 'a' as a temp to save u0, for the score test*/
	a[i] = u[i];

    *flag= cholesky2(imat, nvar, *tol_chol);
    chsolve2(imat,nvar,a);        /* a replaced by  a *inverse(i) */

    *sctest=0;
    for (i=0; i<nvar; i++)
	*sctest +=  u[i]*a[i];

    /*
    **  Never, never complain about convergence on the first step.  That way,
    **  if someone HAS to they can force one iter at a time.
    */
    for (i=0; i<nvar; i++) {
	newbeta[i] = beta[i] + a[i];
	}
    if (*maxiter==0) {
	chinv2(imat,nvar);
	for (i=1; i<nvar; i++)
	    for (j=0; j<i; j++)  imat[i][j] = imat[j][i];
	*flag=0;
	return;   /* and we leave the old beta in peace */
	}

    /*
    ** here is the main loop
    */
    halving =0 ;             /* =1 when in the midst of "step halving" */
    for (iter=1; iter<=*maxiter; iter++) {
	newlk =0;
	for (i=0; i<nvar; i++) {
	    u[i] =0;
	    for (j=0; j<nvar; j++)
		imat[i][j] =0;
	    }

	for (person=0; person<n; person++) {
	    zbeta = 0;      /* form the term beta*z   (vector mult) */
	    for (i=0; i<nvar; i++)
		zbeta += newbeta[i]*covar[i][person];
	    score[person] = exp(zbeta + offset[person]);
	    }

	for (person=0; person<n; ) {
	    R_CheckUserInterrupt();
	    if (event[person]==0) person++;
	    else {
		denom =0;
		for (i=0; i<nvar; i++) {
		    a[i] =0;
		    for (j=0; j<nvar; j++)  cmat[i][j]=0;
		    }

		/*
		** Compute whom is in the risk group, and #deaths
		*/
		nrisk=0;
		deaths=0;
		time = stop[person];
		for (k=person; k<n; k++) {
		    if (stop[k]==time) deaths += event[k];
		    if (start[k] < time) {
			atrisk[nrisk]=k;
			nrisk++;
			}
		    if (strata[k]==1) break;
		    }

		/*
		** compute the mean and covariance over the risk set (a and c)
		**   It's fast if #deaths=1
		*/
		if (deaths==1) {
		    for (l=0; l<nrisk; l++) {
			k = atrisk[l];
			weight = score[k];
			denom += weight;
			for (i=0; i<nvar; i++) {
			    a[i] = a[i] + weight*covar[i][k];
			    for (j=0; j<=i; j++)
				cmat[i][j] += weight*covar[i][k]*covar[j][k];
			    }
			 }
		    }
		else {
		    /*
		    ** for each unique subset of size "deaths" from the risk set,
		    **  the "new variable" is the sum of x over that set.  It's
		    **  weight is the product of the weights.  I want a to contain
		    **  the weighted sum and c the weighted ss of this new var.
		    */
		    init_doloop(0,nrisk);
		    while(doloop(deaths, index) >=0) {
			for (i=0; i<nvar; i++) newvar[i]=0;
			weight =1;
			for (l=0; l<deaths; l++) {
			    k = atrisk[index[l]];
			    weight *= score[k];
			    for (i=0; i<nvar; i++)  newvar[i]+= covar[i][k];
			    }
			denom += weight;
			for (i=0; i<nvar; i++) {
			    a[i] = a[i] + weight*newvar[i];
			    for (j=0; j<=i; j++)
				cmat[i][j] += weight*newvar[i]*newvar[j];
			    }
			 }
		    }
		/*
		** Add results into u and imat
		*/
		newlk -= log(denom);
		for (i=0; i<nvar; i++) {
		    u[i] -=  a[i]/denom;
		    for (j=0; j<=i; j++)
			imat[j][i] += (cmat[i][j] - a[i]*a[j]/denom)/denom;
		    }
		for (k=person; k<n && stop[k]==time; k++) {
		    if (event[k]==1) {
			newlk +=  log(score[k]);
			for (i=0; i<nvar; i++) u[i] += covar[i][k];
			}
		    person++;
		    if (strata[k]==1) break;
		    }
		}
	    }   /* end  of accumulation loop */

	/* am I done?
	**   update the betas and test for convergence
	*/
	*flag = cholesky2(imat, nvar, *tol_chol);

	if (fabs(1-(loglik[1]/newlk))<=*eps && halving==0) { /* all done */
	    loglik[1] = newlk;
	    chinv2(imat, nvar);     /* invert the information matrix */
	    for (i=1; i<nvar; i++)
		for (j=0; j<i; j++)  imat[i][j] = imat[j][i];
	    for (i=0; i<nvar; i++)
		beta[i] = newbeta[i];
	    *maxiter = iter;
	    return;
	    }

	if (iter==*maxiter) break;  /*skip the step halving and etc */

	if (newlk < loglik[1])   {    /*it is not converging ! */
		halving =1;
		for (i=0; i<nvar; i++)
		    newbeta[i] = (newbeta[i] + beta[i]) /2; /*half of old increment */
		}
	    else {
		halving=0;
		loglik[1] = newlk;
		chsolve2(imat,nvar,u);

		j=0;
		for (i=0; i<nvar; i++) {
		    beta[i] = newbeta[i];
		    newbeta[i] = newbeta[i] +  u[i];
		    }
		}
	}   /* return for another iteration */

    loglik[1] = newlk;
    chinv2(imat, nvar);
    for (i=1; i<nvar; i++)
	for (j=0; j<i; j++)  imat[i][j] = imat[j][i];
    for (i=0; i<nvar; i++)
	beta[i] = newbeta[i];
    *flag= 1000;
    }
