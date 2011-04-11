/*
** Anderson-Gill formulation of the Cox Model, using smart subsets
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
**       nstrat       :number of strata
**       strata(nstrat):sizes of the strata
**       sort1        :sort order for the obs, using stop time, last to first
**       sort2        :sort order for the obs, using start time, last to first
**       offset(n)    :linear offset
**       weights(n)   :case weights
**       eps          :tolerance for convergence.  Iteration continues until
**                       the percent change in loglikelihood is <= eps.
**       tol_chol     : tolerance for the Cholesky routine
**
**  returned parameters
**       means(nv)    :column means of the X matrix
**       beta(nv)     :the vector of answers (at start contains initial est)
**       u(nv)        :score vector
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
**       a2(nvar)
**       cmat(nvar,nvar)       ragged array
**       cmat2(nvar,nvar)
**       newbeta(nvar)         always contains the "next iteration"
**
**  the 6 arrays score, a, cmat, and newbeta are passed as a single
**    vector of storage, and then broken out.
**
**  calls functions:  cholesky2, chsolve2, chinv2
*/
#include <math.h>
#include "survS.h" 
#include "survproto.h"

void agfit3( Sint   *maxiter,  Sint   *nusedx,  Sint   *nvarx, 
	     double *start,    double *stop,    Sint   *event, 
	     double *covar2,   double *offset,  double *weights,
	     Sint   *nstrat,   Sint   *strata,  Sint   *sort1,
	     Sint   *sort2,    double *means,   double *beta, 
	     double *u,        double *imat2,   double loglik[2], 
	     Sint   *flag,     double *work,   
	     double *eps,      double *tol_chol, double *sctest)
{

    int i,j,k,person;
    int indx2, istrat, p;
    int ksave;
    int     iter;
    int     nused, nvar;

    double **covar, **cmat, **imat;  /*ragged array versions*/
    double *a, *newbeta;
    double *a2, **cmat2;
    double *score;
    double  denom, zbeta, risk;
    double  time;
    double  temp, temp2;
    double  newlk =0;
    int     halving;    /*are we doing step halving at the moment? */
    double     method;
    double  meanwt;
    int itemp, deaths;
    double efron_wt, d2;

    nused = *nusedx;
    nvar  = *nvarx;
    method= *sctest;
    /*
    **  Set up the ragged arrays
    */
    covar= dmatrix(covar2, nused, nvar);
    imat = dmatrix(imat2,  nvar, nvar);
    cmat = dmatrix(work,   nvar, nvar);
    cmat2= dmatrix(work + nvar*nvar, nvar, nvar);
    a = work + 2*nvar*nvar;
    a2= a+nvar;
    newbeta = a2 + nvar;
    score   = newbeta + nvar;

    /*
    ** Subtract the mean from each covar, as this makes the regression
    **  much more stable
    */
    for (i=0; i<nvar; i++) {
	temp=0;
	for (person=0; person<nused; person++) temp += covar[i][person];
	temp /= nused;
	means[i] = temp;
	for (person=0; person<nused; person++) covar[i][person] -=temp;
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

    for (person=0; person<nused; person++) {
	zbeta = 0;      /* form the term beta*z   (vector mult) */
	for (i=0; i<nvar; i++)
	    zbeta += beta[i]*covar[i][person];
	score[person] = coxsafe(zbeta + offset[person]);
        }

    /*
    **  'person' walks through the stop times from largest to smallest
    **     (sort1[0] points to the largest stop time, sort1[1] the next, ...)
    **  'time' is the time of current interest, which also goes from smallest
    **      to largest.  
    **  'indx2' walks through the start times.  It will be smaller than 
    **    'person': if person=27 that means that 27 subjects have stop >=time,
    **    and are thus potential members of the risk set.  If 'indx2' =9,
    **    that means that 9 subjects have start >=time and thus are NOT part
    **    of the risk set.  (stop > start for each subject guarrantees that
    **    the 9 are a subset of the 27). 
    **  Basic algorithm: move 'person' forward, adding the new subject into
    **    the risk set.  If this is a new, unique death time, take selected
    **    old obs out of the sums, add in obs tied at this time, then
    **    add terms to the loglik, etc.
    */
    istrat=0;
    indx2 =0;
    denom =0;
    for (i=0; i<nvar; i++) {
	a[i] =0;
	for (j=0; j<nvar; j++) {
	    cmat[i][j]=0;
	    }
	}
    
    for (person=0; person<nused;) {
	p = sort1[person];
	if (event[p]==0){
	    risk = exp(score[p]) * weights[p];
	    denom += risk;
	    for (i=0; i<nvar; i++) {
		a[i] += risk*covar[i][p];
		for (j=0; j<=i; j++)
		    cmat[i][j] += risk*covar[i][p]*covar[j][p];
		}
	    person++;
	    }
	else {
	    time = stop[p];
	    /*
	    ** subtract out the subjects whose start time is to the right
	    */
	    for (; indx2<strata[istrat]; indx2++) {
		p = sort2[indx2];
		if (start[p] < time) break;
		risk = exp(score[p]) * weights[p];
		denom -= risk;
		for (i=0; i<nvar; i++) {
		    a[i] -= risk*covar[i][p];
		    for (j=0; j<=i; j++)
			cmat[i][j] -= risk*covar[i][p]*covar[j][p];
		    }
		}
	    /*
	    ** compute the averages over subjects with
	    **   exactly this death time (a2 & c2)
	    ** (and add them into a and cmat while we're at it)
	    */
	    efron_wt =0;
	    meanwt =0;
	    for (i=0; i<nvar; i++) {
		a2[i]=0;
		for (j=0; j<nvar; j++) {
		    cmat2[i][j]=0;
		    }
		}
	    deaths=0;
	    for (k=person; k<strata[istrat]; k++) {
		p = sort1[k];
		if (stop[p] < time) break;
		risk = exp(score[p]) * weights[p];
		denom += risk;
		for (i=0; i<nvar; i++) {
		    a[i] += risk*covar[i][p];
		    for (j=0; j<=i; j++)
			cmat[i][j] += risk*covar[i][p]*covar[j][p];
		    }
		if (event[p]==1) {
		    deaths += event[p];
		    efron_wt += risk*event[p];
		    meanwt += weights[p];
		    for (i=0; i<nvar; i++) {
			a2[i]+= risk*covar[i][p];
			for (j=0; j<=i; j++)
			    cmat2[i][j] += risk*covar[i][p]*covar[j][p];
			}
		     }
		}
	    ksave = k;
	    /*
	    ** Add results into u and imat for all events at this time point
	    */
	    meanwt /= deaths;
	    itemp = -1;
	    for (; person<ksave; person++) {
		p = sort1[person];
		if (event[p]==1) {
		    itemp++;
		    temp = itemp*method/deaths;
		    d2 = denom - temp*efron_wt;
		    loglik[1] +=  weights[p]*score[p] -meanwt *log(d2);
		    for (i=0; i<nvar; i++) {
			temp2 = (a[i] - temp*a2[i])/d2;
			u[i] += weights[p]*covar[i][p] - meanwt*temp2;
			for (j=0; j<=i; j++)
			    imat[j][i] += meanwt* (
					(cmat[i][j] - temp*cmat2[i][j])/d2-
					   temp2*(a[j]-temp*a2[j])/d2);
			}
		    }
		}
	    }

	if (person == strata[istrat]) {
	    istrat++;
	    denom =0;
	    indx2 = person;
	    for (i=0; i<nvar; i++) {
		a[i] =0;
		for (j=0; j<nvar; j++) {
		    cmat[i][j]=0;
		    }
		}
	    }
	}   /* end  of accumulation loop */

    loglik[0] = loglik[1];   /* save the loglik for iteration zero  */

    /* am I done?
    **   update the betas and test for convergence
    */

    for (i=0; i<nvar; i++) /*use 'a' as a temp to save u0, for the score test*/
	a[i] = u[i];

    *flag = cholesky2(imat, nvar, *tol_chol);
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

	for (person=0; person<nused; person++) {
	    zbeta = 0;      /* form the term beta*z   (vector mult) */
	    for (i=0; i<nvar; i++)
		zbeta += newbeta[i]*covar[i][person];
	    score[person] = coxsafe(zbeta + offset[person]);
#if 0
	    /* I belive this text is no longer needed due to the coxsafe call
	     */
	    if (zbeta > 20 && *maxiter>1) {
		/*
		** If the above happens, then 
		**   1. There is a real chance for catastrophic cancellation
		**       in the computation of "denom", which leads to
		**       numeric failure via log(neg number) -> inf loglik
		**   2. A risk score for one person of exp(20) > 400 million
		**       is either an infinite beta, in which case any
		**       reasonable coefficient will do, or a big overreach
		**       in the Newton-Raphson step.
		** In either case, a good solution is step halving.  However,
		**   if the user asked for exactly 1 iteration, return it.
		** 
		** Why 20?  Most machines have about 16 digits of precision,
		**   and this preserves approx 7 digits in the subtraction
		**   when a high risk score person leaves the risk set.
		**   (Because of centering, the average risk score is about 0).
		**   Second, if eps is small and beta is infinite, we rarely
		**   get a value above 16.  So a 20 is usually a NR overreach.
		** A data set with zbeta=54 on iter 1 led to this fix, the
		**   true final solution had max values of 4.47.    
		*/
		halving=1;
		for (i=0; i<nvar; i++)
		    newbeta[i] = (newbeta[i] + beta[i])/2.1;
		person = -1;  /* force the loop to start over */
		}
#endif
	    }

	istrat=0;
	indx2 =0;
	denom =0;
	for (i=0; i<nvar; i++) {
	    a[i] =0;
	    for (j=0; j<nvar; j++) {
		cmat[i][j]=0;
		}
	    }
    
	for (person=0; person<nused;) {
	    p = sort1[person];
	    if (event[p]==0){
		risk = exp(score[p]) * weights[p];
		denom += risk;
		for (i=0; i<nvar; i++) {
		    a[i] += risk*covar[i][p];
		    for (j=0; j<=i; j++)
			cmat[i][j] += risk*covar[i][p]*covar[j][p];
		    }
		person++;
		}
	    else {
		time = stop[p];
		/*
		** subtract out the subjects whose start time is to the right
		*/
		for (; indx2<strata[istrat]; indx2++) {
		    p = sort2[indx2];
		    if (start[p] < time) break;
		    risk = exp(score[p]) * weights[p];
		    denom -= risk;
		    for (i=0; i<nvar; i++) {
			a[i] -= risk*covar[i][p];
			for (j=0; j<=i; j++)
			    cmat[i][j] -= risk*covar[i][p]*covar[j][p];
			}
		    }

		/*
		** compute the averages over this death time (a2 & c2)
		*/
		efron_wt =0;
		meanwt =0;
		for (i=0; i<nvar; i++) {
		    a2[i]=0;
		    for (j=0; j<nvar; j++) {
			cmat2[i][j]=0;
			}
		    }
		deaths=0;
		for (k=person; k<strata[istrat]; k++) {
		    p = sort1[k];
		    if (stop[p] < time) break;
		    risk = exp(score[p]) * weights[p];
		    denom += risk;
		    for (i=0; i<nvar; i++) {
			a[i] += risk*covar[i][p];
			for (j=0; j<=i; j++)
			    cmat[i][j] += risk*covar[i][p]*covar[j][p];
			}
		    if (event[p]==1) {
			deaths += event[p];
			efron_wt += risk*event[p];
			meanwt += weights[p];
			for (i=0; i<nvar; i++) {
			    a2[i]+= risk*covar[i][p];
			    for (j=0; j<=i; j++)
				cmat2[i][j] += risk*covar[i][p]*covar[j][p];
			    }
			}
		    }
		ksave = k;

		/*
		** Add results into u and imat 
		*/
		meanwt /= deaths;
		itemp = -1;
		for (; person<ksave; person++) {
		    p = sort1[person];
		    if (event[p]==1) {
			itemp++;
			temp = itemp*method/deaths;
			d2 = denom - temp*efron_wt;
			newlk +=  weights[p]*score[p] -meanwt *log(d2);

			for (i=0; i<nvar; i++) {
			    temp2 = (a[i] - temp*a2[i])/d2;
			    u[i] += weights[p]*covar[i][p] - meanwt*temp2;
			    for (j=0; j<=i; j++)
				imat[j][i] += meanwt* (
				    (cmat[i][j] - temp*cmat2[i][j])/d2-
				    temp2*(a[j]-temp*a2[j])/d2);
			    }
			}
		    }
		}
	
	    if (person == strata[istrat]) {
		istrat++;
		denom =0;
		indx2 = person;
		for (i=0; i<nvar; i++) {
		    a[i] =0;
		    for (j=0; j<nvar; j++) {
			cmat[i][j]=0;
			}
		    }
		}
	    }   /* end  of accumulation loop */

	/* am I done?
	**   update the betas and test for convergence
	*/
	*flag = cholesky2(imat, nvar, *tol_chol);

	if (fabs(1-(loglik[1]/newlk))<=*eps  && halving==0) { /* all done */
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
