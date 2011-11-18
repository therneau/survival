/*
** Create the cohort survival curve(s) for a set of subjects.
**
**   This is similar in output to pyears3.  However, the hard part now is
**        computation of a variance.
**
**  Input
**    n=# of subjects
**    nvar - number of vars in xdata
**    ncurve - number of curves that will be produced
**    npt    - number of points on each curve, e.g, #unique death times
**    se     - 0:no standard errors  1:do standard errors
**    score  - the risk score for the old subjects
**    y -  the max fu time for each subject
**    score[n] - vector of weights
**    strata = the group designator (for multiple survival curves)
**    r  =    an nvar column matrix containing the X data
**              for the new subjects.
**    var    = Cox variance matrix
**    xmean  = vector of means of the data, for centering.
**    cn     =  n for the original (Cox) data
**    cy     =  3 column matrix of original data
**    cx     =  the original X data -- will be a dummy is se =0
**    method   = the methods
**                    individual:  1=Kalbfleisch/Prentice  2= Tsiatis/Breslow
**                              3= Tsiatis, Efron approx
**                    cohort: + (10 if conditional)  (0 otherwise)
**
** Output
**    surv[ncurve][npt]  the conditional survival at t
**    varh[ncurve][npt]  variance of the survival (if requested)
**    used[ncurve][npt]  # of subjects contributing to each time point
**
**  Work space is allocated as needed
**    for the calculations on the "old" data
**          a[nvar], a2[nvar]
**          score[cn]  -- the risk score, per subject
**    and on the "new"
**          newscore[n] -- risk score
**          tvar[n][n]  -- holds cross product terms -- only if nvar>0
**          isurv[n]    -- individual survival curves
**
**  cx and cy must be sorted by (event before censor) within stop time
*/
#include <math.h>
#include "survS.h"
#include "survproto.h"

static double   *y,
		*nscore,
		**newx,
		**surv,
		**vsurv,
		*isurv,
		**used,
		**tvar;
static int      *strata;
static double   ttime,    /* Some HP compilers choke on "time" as a variable */
		**imat,
		*mean;
static int      death,
		ncurve,
		se,
		nvar,
		n;
static void    addup();

void agsurv3(Sint   *sn,    Sint   *snvar,    Sint   *sncurve, 
	     Sint   *snpt,  Sint   *sse,      double *score, 
	     double *sy,    Sint   *grpx,     double *r,        double *coef, 
	     double *var,   double *xmean,    Sint   *scn, 
	     double *cy,    double *cx,       double *ssurv,
	     double *varh,  double *sused,    Sint   *smethod)
{

    int i,j,k,l;
    double *start, *stop, *event;
    int cn;
    int npt,
	nvar2,
	method;
    int kk=0, psave;
    int itime;
    int person;
    int deaths, nrisk;
    int need;
    double *a=0, *a2=0;
    double weight=0,
	   e_denom,
	   denom;
    double inc,
	   sumt,
	   km =0;
    double temp,
	   downwt,
	   d2;
    double haz,
	   varhaz;
    double **oldx =0;


    n = *sn;  nvar = *snvar;
    cn = *scn; npt = *snpt;
    se = *sse;
    ncurve = *sncurve;
    method = *smethod;
    death = method/10;
    method = method - death*10;
    y = sy;
    start = cy;
    stop  = cy+ cn;
    event = cy+ cn+ cn;
    strata = grpx;

    /*
    ** scratch space
    */
    need = 2*n + se*nvar*(2+ n*(n+1)/2) + nvar;
    nscore = (double *) ALLOC(need, sizeof(double));
    for (i=0; i<need; i++) nscore[i] =0.0;  /* R doesn't zero the memory */
    mean = nscore +n;
    isurv  = mean + nvar;
    for (i=0; i<n; i++) isurv[i]=1;
    if (se==1) {
	a = isurv + n;
	a2= a + nvar;
	tvar = (double **) ALLOC(n, sizeof(double *));
	/* be tricky here, as only the bottom half is used */
	tvar[0] = a2 + nvar;
	for (i=1; i<n; i++)
	    tvar[i] = tvar[i-1] +i;
	}

    /*
    **  Set up the ragged arrays
    */
    if (se==1) oldx = dmatrix(cx, cn, nvar);
    newx = dmatrix(r, n, nvar);
    imat = dmatrix(var,  nvar, nvar);
    surv = dmatrix(ssurv, npt, ncurve);
    vsurv = dmatrix(varh,  npt, ncurve);
    used = dmatrix(sused, npt, ncurve);

    for (i=0; i<ncurve; i++)
	for (j=0; j<npt; j++)  surv[i][j] =1;

    /*
    ** compute the risk scores
    */
    for (i=0; i<n; i++) {
	nscore[i] =0;
	for (j=0; j<nvar; j++) {
	    nscore[i] += coef[j]* (newx[j][i] - xmean[j]);
	    }
	nscore[i] = exp(nscore[i]);
	}

    /*
    ** outer loop keeps a running track of the baseline hazard
    **  the addup function adds things up
    */
    itime =0;
    nvar2 = nvar * se;    /*simpler than a lot of " if (se==1)" statements */
    for (person=0; person<cn;) {
	if (event[person]==0) person++;
	else {
	    /*
	    ** compute the mean and denominator over the risk set
	    */
	    denom =0;
	    e_denom=0;
	    for(i=0; i<nvar2; i++){
		a[i] =0;
		a2[i]=0;
		}
	    ttime = stop[person];
	    nrisk =0;
	    deaths=0;
	    for (k=person; k<cn; k++) {
		if (start[k] < ttime) {
		    nrisk++;
		    weight = score[k];
		    denom += weight;
		    for (i=0; i<nvar2; i++) {
			a[i] += weight*(oldx[i][k] - xmean[i]);
			}
		     }
		if (stop[k]==ttime && event[k]==1) {
		    kk=k;
		    deaths++;
		    e_denom += weight;
		    for (i=0; i<nvar2; i++) {
			a2[i] += weight*(oldx[i][k] - xmean[i]);
			}
		    }
		}

	    /*
	    ** Now compute the increment in the hazard and variance at "time"
	    */
	    if (method <3) for (i=0; i<nvar2; i++) mean[i] = a[i]/denom;
	    if (method==1) {
		for (psave=person; psave<cn && stop[psave]==ttime; psave++) 
		/*
		** kalbfleisch estimator requires iteration;
		*/
		if (deaths == nrisk) km=0;
		else if (deaths ==1) {
		    km = pow(1- score[kk]/denom, 1/score[kk]);
		    }
		else {           /*find the zero of an equation */
		    km = .5;
		    inc = .25;
		    for (l=0; l<35; l++) { /* bisect it to death */
			sumt =0;
			for (k=person; k<psave; k++) {
			    if (event[k] ==1)
				sumt +=  score[k]/(1-pow(km, score[k]));
			    }
			if (sumt < denom)  km += inc;
			     else          km -= inc;
			inc = inc/2;
			}
		    }
		 if (km==0)
		    addup(itime, 0.0, 0.0);
		 else {
		    haz = log(km);
		    varhaz = deaths/(denom *(denom-e_denom));  /* Greenwood */
		    addup(itime, haz, varhaz);
		    }
		 person = psave;
		 }

	    else if (method==2) {
		addup(itime, deaths/denom, deaths/(denom*denom));
		for (; person<cn && stop[person]==ttime; person++);
		}
	    else {
		temp =0;  haz=0; varhaz=0;
		for (k=person; k<cn && stop[k]==ttime; k++) {
		    if (event[k]==1) {
			downwt = temp++/deaths;
			d2 = (denom - downwt*e_denom);
			haz = 1/d2;
			varhaz = 1/(d2*d2);
			for (i=0; i<nvar2; i++)
			    mean[i] = (a[i] - downwt*a2[i])/ d2;
			addup(itime, haz, varhaz);
			}
		    person++;
		    }
		}
	    itime++;
	    }
	}
    }


static void addup(itime, haz, var)
int itime;
double haz, var;
    {
    int i, j, k, l;
    int     pstart,
	    ic;
    double  temp,
	    totsurv,
	    totvar,
	    wt,
	    nn;

    if (var==0) {
	/*
	** KM method, and everyone in the referent group died
	*/
	for (i=0; i<ncurve; i++) {
	    surv[i][itime] =0;
	    if (nvar>0) vsurv[i][itime]=0;
	    }
	return;
	}

    /*
    ** Note that the subjects are sorted in strata order
    */
    pstart=0;
    for (ic=0; ic<ncurve; ic++) {
	nn =0;
	wt =0;
	totsurv=0;
	totvar =0;
	for (i=pstart; i<n && strata[i]==ic; i++) {
	    nn++;
	    if (y[i] >= ttime) {
		temp =  -haz*nscore[i];  /*increment to the individual hazard*/
		if  (death==0) {
		    wt += isurv[i];
		    totsurv += exp(temp) * isurv[i];
		    }
		else {
		    wt += 1;
		    totsurv += temp;
		    }
		isurv[i] *= exp(temp);
		}
	    /*
	    ** The variance is computed as though it were the Ederer est, always
	    */
	    if (se==1)  {  /* Do the variance term (nasty) */
		for (j=pstart; j<=i; j++) {
		    temp =0;
		    for (k=0; k<nvar; k++) {
			temp += (newx[k][i]-mean[k])*(newx[k][j]-mean[k])*
					  imat[k][k];
			for (l=0; l<k; l++)
			    temp += ((newx[k][i]-mean[k])*(newx[l][j]-mean[l]) +
				     (newx[k][j]-mean[k])*(newx[l][i]-mean[l]))*
				     imat[k][l];
			 }
		    tvar[i][j] += (1+ temp) * var;
		    temp = nscore[i]*nscore[j]* tvar[i][j] *
				 isurv[i] * isurv[j];
		    if (i==j) totvar += temp;
		    else      totvar += temp + temp;
		    }
		}
	    }
	/*
	** save the answer
	**  surv is "accumulated" because with the Efron estimate multiple
	**  calls are made to addup() with the same value of itime.
	*/
	used[ic][itime] = nn;
	if (death==0) surv[ic][itime] *= totsurv/wt;
	else          surv[ic][itime] *= exp(totsurv/wt);
	if (se==1) vsurv[ic][itime] = totvar/(nn*nn);
	pstart =i;
	}
    }
