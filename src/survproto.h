/*
** $Id $
**
** Prototypes of all the survival functions
**  Including this in each routine helps prevent mismatched argument errors
*/

void agexact(Sint *maxiter,  Sint *nusedx,   Sint *nvarx,   double *start, 
	     double *stop,   Sint *event,    double *covar2,double *offset, 
	     Sint   *strata, double *means,  double *beta,  double *u, 
	     double *imat2,  double loglik[2], Sint *flag,  double *work, 
	     Sint   *work2,  double *eps,    double *tol_chol, double *sctest);

void agfit2( Sint   *maxiter,  Sint   *nusedx,  Sint   *nvarx, 
	     double *start,    double *stop,    Sint   *event, 
	     double *covar2,   double *offset,  double *weights,
	     Sint   *strata,   double *means,   double *beta, 
	     double *u,        double *imat2,   double loglik[2], 
	     Sint   *flag,     double *work,    Sint   *end,
	     double *eps,      double *tol_chol,double *sctest);

void agfit4_a(Sint *nusedx, Sint *nvarx, double *yy, 
	       double *covar2, double *offset2,
	       double *weights2, Sint *strata2,
	       double *means, double *beta, double *u, 
	       double *loglik, 
	       Sint *methodx, Sint *ptype2, Sint *pdiag2,
	       Sint *nfrail,  Sint *frail2);

void agfit4_b(Sint *maxiter, Sint *nusedx, Sint *nvarx, 
	       double *beta, double *u,
	       double *imat2,  double *jmat2, double *loglik, 
	       Sint *flag,  double *eps, double *tolerch, Sint *methodx, 
	       Sint *nfrail, double *fbeta, double *fdiag);

void agfit_null(Sint   *n,      Sint   *method,   double *start, double *stop, 
		Sint   *event,  double * offset,  double *weights,
		Sint   *strata, double loglik[2]);

void aghaz2(Sint   *n,     double *start,   double *stop,   Sint   *event, 
	    double *score, Sint   * strata, double *hazard, double * cumhaz);

void agmart(Sint   *n,     Sint   *method,  double *start,   double *stop, 
	    Sint   *event, double *score,   double *wt,      Sint   *strata, 
	    double *resid);

void agres12(Sint   *nx,     Sint   *nvarx,   double *y,    double *covar2, 
	     Sint   *strata, double *score,   Sint *method, double *resid2, 
	     double *a);

void agscore(Sint   *nx,       Sint   *nvarx,      double *y,
	     double *covar2,   Sint   *strata,     double *score,
	     double *weights,  Sint   *method,     double *resid2, double *a);

void agsurv1(Sint   *sn,     Sint   *snvar,  double *y,      double *score, 
	     Sint   *strata, 
	     double *surv,   double *varh,   Sint   *snsurv,
	     double *xmat,   double *d,      double *varcov, double *yy,
	     Sint   *shisn,  double *hisy,   double *hisxmat,double *hisrisk, 
	     Sint   *hisstrat);

void agsurv2(Sint   *sn,      Sint   *snvar,    double *y, 
	     double *score,   Sint   *strata,   double *wt,    double *surv, 
	     double *varh,    double *xmat,     double *varcov, 
	     Sint   *snsurv,  double *d,        Sint   *sncurve,
             double *newx,    double *newrisk);

void agsurv3(Sint   *sn,    Sint   *snvar,    Sint   *sncurve, 
	     Sint   *snpt,  Sint   *sse,      double *score, 
	     double *sy,    double *r,        double *coef, 
	     double *var,   double *cmean,    Sint   *scn, 
	     double *cy,    double *cx,       double *ssurv,
	     double *varh,  double *sused,    Sint   *smethod);

void chinv2  (double **matrix, int n);
int cholesky2(double **matrix, int n, double toler);
void chsolve2(double **matrix, int n, double *y);
void chinv3(double **matrix , int n, int m, double *fdiag);
int cholesky3(double **matrix, int n, int m, double *diag, double toler);
void chsolve3(double **matrix, int n, int m, double *diag, double *y);

void coxdetail(Sint   *nusedx,   Sint   *nvarx,    Sint   *ndeadx, 
	       double *y,        double *covar2,   Sint   *strata,  
	       double *score,    double *weights,  double *means2, 
	       double *u2,       double *var,      Sint   *rmat,
	       double *nrisk2,   double *work);

void coxfit2(Sint   *maxiter,   Sint   *nusedx,    Sint   *nvarx, 
	     double *time,      Sint   *status,    double *covar2, 
	     double *offset,	double *weights,   Sint   *strata,
	     double *means,     double *beta,      double *u, 
	     double *imat2,     double loglik[2],  Sint   *flag, 
	     double *work,	double *eps,       double *tol_chol,
	     double *sctest);

void coxfit4_a(Sint *nusedx, Sint *nvarx, double *yy, 
               double *covar2, double *offset2,
               double *weights2, Sint *strata2,
               double *means, double *beta, double *u, 
               double *loglik, 
               Sint *methodx, Sint *ptype2, Sint *pdiag2,
               Sint *nfrail,  Sint *frail2);

void coxfit4_b(Sint *maxiter, Sint *nusedx, Sint *nvarx, 
               double *beta, double *u,
               double *imat2,  double *jmat2, double *loglik, 
               Sint *flag,  double *eps, double *tolerch, Sint *methodx, 
               Sint *nfrail, double *fbeta, double *fdiag);

void coxfit4_c (Sint *nusedx, Sint *nvar, Sint *methodx, double *expect);

void coxfit_null(Sint   *nusedx,    Sint   *method,   double *time, 
		 Sint   *status,    double *score,    double *weights, 
		 Sint   *strata,    double *loglik, double *resid);

void coxhaz2(Sint   *n,      double *score,   Sint   *mark, 
	     Sint   *strata, double *hazard,  double *cumhaz);

void coxmart(Sint   *sn,     Sint   *method,    double *time, 
	     Sint   *status, Sint   * strata,   double *score, 
	     double *wt,     double *expect);

void coxph_wtest(Sint *nvar2, Sint *ntest, double *var, double *b,
                 double *scratch, double *tolerch);

void coxscho(Sint   *nusedx,    Sint   *nvarx,    double *y, 
	     double *covar2,    double *score,    Sint   *strata,  
	     Sint   *method2,   double *work);

void coxscore(Sint   *nx,      Sint   *nvarx,    double *y, 
	      double *covar2,  Sint   *strata,   double *score, 
	      double *weights, Sint   *method,   double *resid2,
	      double *scratch);

double **dmatrix(double *array, int ncol, int nrow);

void init_doloop(int min, int max);
int doloop      (int nloops, int *index);

void pyears1(Sint   *sn,      Sint   *sny,      Sint   *sdoevent, 
	     double *sy,      double *wt,       
	     Sint   *sedim,   Sint   *efac, 
	     Sint   *edims,   double *secut,    double *expect, 
	     double *sedata,  Sint   *sodim,    Sint   *ofac, 
	     Sint   *odims,   double *socut,    Sint   *smethod, 
	     double *sodata,  double *pyears,   double *pn, 
	     double *pcount,  double *pexpect,  double *offtable);

void pyears2(Sint   *sn,      Sint   *sny,   Sint   *sdoevent, 
	     double *sy,      double *wt,    
	     Sint   *sodim,   Sint   *ofac, 
	     Sint   *odims,   double *socut, double *sodata,
	     double *pyears,  double *pn,    double *pcount, 
	     double *offtable);

void pyears3(Sint   *sdeath,    Sint   *sn,    Sint   *sedim, 
	     Sint   *efac,      Sint   *edims, double *secut, 
	     double *expect,    double *sx,    double *y, 
	     Sint   *sntime,    Sint   *sngrp, double *times,
	     double *esurv,     Sint   *nsurv);

double pystep(int nc,        int  *index,  int  *index2,   double *wt, 
	      double *data,  Sint *fac,    Sint *dims,     double **cuts, 
	      double step,   int  edge);

void survdiff2(Sint   *nn,     Sint   *nngroup,    Sint   *nstrat, 
	       double *rho,    double *time,       Sint   *status, 
	       Sint   *group,  Sint   *strata,	   double *obs, 
	       double *exp,    double *var,        double *risk, 
	       double *kaplan);

void survfit2(Sint   *sn,      double *y,       double *wt,
	      Sint   *strata,  Sint   *method,  Sint   *error, 
	      double *mark,    double *surv,	double *varh,
	      double *risksum);

void survfit3(Sint   *sn,        double *y,               double *wt,
	      Sint   *strata,    Sint   *method,          Sint   *error, 
	      Sint   *nstrat,    double *ntimes_strata,  
	      double *timelist,  double *weighted_event,  double *surv,
	      double *varh,	 double *risksum,         double *enter,
	      double *exit_censored);


SEXP survreg6(SEXP maxiter2,   SEXP nvarx,  SEXP y,
	      SEXP ny2,        SEXP covar2, SEXP wtx,
	      SEXP offset2,    SEXP beta2,  SEXP nstratx,
	      SEXP stratax,    SEXP epsx,   SEXP tolx,       
	      SEXP dist,       SEXP expr,   SEXP rho);

double survregc1(int n,          int nvar,     int nstrat,      int whichcase,
		 double *beta,   int dist,     Sint *strat,     double *offset,
		 double *time1,  double *time2, double *status, double *wt,
		 double **covar, double **imat, double **JJ,    double *u, 
		 SEXP expr,      SEXP rho,      double *dummy,  int nf,
		 Sint *frail,    double *fdiag, double *jdiag );

double survregc2(int n,          int nvar,     int nstrat,      int whichcase,
		 double *beta,   int dist,     Sint *strat,     double *offset,
		 double *time1,  double *time2, double *status, double *wt,
		 double **covar, double **imat, double **JJ,    double *u, 
		 SEXP expr,      SEXP rho,      double *dummy,  int nf,
		 Sint *frail,    double *fdiag, double *jdiag );

void survpenal(int whichcase, int nfrail,    int  nvar2,    double **hmat, 
	       double **JJ,   double *hdiag, double *jdiag,
	       double *u,     double *beta,  double *loglik,
	       int ptype,     int pdiag,     SEXP pexpr1,   double *cptr1, 
	       SEXP pexpr2,   double *cptr2, SEXP rho);

void coxpenal (int whichcase, int nfrail,   int  nvar,     double **hmat, 
	       double *hdiag, double *u,    double *fbeta, double *beta,  
	       double *penalty,
	       int ptype,     int pdiag,     SEXP pexpr1,   double *cptr1, 
	       SEXP pexpr2,   double *cptr2, SEXP rho);
