/*
** Prototypes of all the survival functions
**  Including this in each routine helps prevent mismatched argument errors
*/
void agfit3( Sint   *maxiter,  Sint   *nusedx,  Sint   *nvarx, 
	     double *start,    double *stop,    Sint   *event, 
	     double *covar2,   double *offset,  double *weights,
	     Sint   *nstrat,   Sint   *strata,  Sint   *sort1,
	     Sint   *sort2,    double *means,   double *beta, 
	     double *u,        double *imat2,   double loglik[2], 
	     Sint   *flag,     double *work,   
	     double *eps,      double *tol_chol, double *sctest);

void agfit5_a(Sint *nusedx,     Sint *nvarx,     double *yy, 
	      double *covar2,   double *offset2, double *weights2, 
	      int   *strata,    Sint   *sort,    double *means,   
              double *beta,     double *u,       double *loglik, 
	      Sint *methodx,    Sint *ptype2,    Sint *pdiag2,
	      Sint *nfrail,     Sint *frail2,
              void *fexpr1,     void *fexpr2,    void *rho) ;

void agfit5_b( Sint *maxiter,   Sint *nusedx,    Sint *nvarx, 
	       Sint *strata,    double *beta,    double *u,
	       double *imat2,   double *jmat2,   double *loglik, 
	       Sint *flag,      double *eps,     double *tolerch, 
	       Sint *methodx,   Sint *nfrail,     double *fbeta, 
	       double *fdiag,
               void *fexpr1,    void *fexpr2,     void *rho);

void agfit5_c(Sint *nusedx,   Sint *nvar,    Sint *strata,
	      Sint *methodx,  double *expect) ;
 
void agexact(Sint *maxiter,  Sint *nusedx,   Sint *nvarx,   double *start, 
	     double *stop,   Sint *event,    double *covar2,double *offset, 
	     Sint   *strata, double *means,  double *beta,  double *u, 
	     double *imat2,  double loglik[2], Sint *flag,  double *work, 
	     Sint   *work2,  double *eps,    double *tol_chol, double *sctest);

void agmart(Sint   *n,     Sint   *method,  double *start,   double *stop, 
	    Sint   *event, double *score,   double *wt,      Sint   *strata, 
	    double *resid);

void agmart2(Sint   *n,     Sint   *method,  double *start,   double *stop, 
	     Sint   *event,  Sint   *nstrat,  Sint *strata,    Sint *sort1,
	     Sint   *sort2,  double *score,   double *wt,      
	     double *resid,  double *haz);
void agscore(Sint   *nx,       Sint   *nvarx,      double *y,
	     double *covar2,   Sint   *strata,     double *score,
	     double *weights,  Sint   *method,     double *resid2, double *a);

void agsurv3(Sint   *sn,    Sint   *snvar,    Sint   *sncurve, 
	     Sint   *snpt,  Sint   *sse,      double *score, 
	     double *sy,    Sint   *grpx,     double *r,        double *coef, 
	     double *var,   double *xmean,    Sint   *scn, 
	     double *cy,    double *cx,       double *ssurv,
	     double *varh,  double *sused,    Sint   *smethod);

void agsurv4(Sint   *ndeath,   double *risk,    double *wt,
             Sint   *sn,        double *denom,   double *km);

void agsurv5(Sint *n2,     Sint *nvar2,  Sint *dd, double *x1,  
             double *x2,   double *xsum, double *xsum2, 
             double *sum1, double *sum2, double *xbar) ;

void chinv2  (double **matrix, int n);
int cholesky2(double **matrix, int n, double toler);
void chsolve2(double **matrix, int n, double *y);
void chinv3(double **matrix , int n, int m, double *fdiag);
int cholesky3(double **matrix, int n, int m, double *diag, double toler);
void chsolve3(double **matrix, int n, int m, double *diag, double *y);

SEXP concordance1(SEXP y, SEXP wt2,  SEXP indx2, SEXP ntree2);
 
SEXP concordance2(SEXP y,     SEXP wt2,  SEXP indx2, SEXP ntree2,
                  SEXP sortstop, SEXP sortstart) ;

SEXP coxcount1(SEXP y2, SEXP strat2) ;
SEXP coxcount2(SEXP y2, SEXP isort1, SEXP isort2, SEXP strat2) ;

void coxdetail(Sint   *nusedx,   Sint   *nvarx,    Sint   *ndeadx, 
	       double *y,        double *covar2,   Sint   *strata,  
	       double *score,    double *weights,  double *means2, 
	       double *u2,       double *var,      Sint   *rmat,
	       double *nrisk2,   double *work);
 
SEXP coxexact(SEXP maxiter2,  SEXP y2, 
              SEXP covar2,    SEXP offset2, SEXP strata2,
              SEXP ibeta,     SEXP eps2,    SEXP toler2) ;
 
void coxfit5_a(Sint *nusedx,     Sint *nvarx,     double *yy, 
 	       double *covar2,   double *offset2, double *weights2, 
	       int   *strata,    Sint   *sort,    double *means,   
               double *beta,     double *u,       double *loglik, 
	       Sint *methodx,    Sint *ptype2,    Sint *pdiag2,
	       Sint *nfrail,     Sint *frail2,
               void *fexpr1,     void *fexpr2,    void *rho) ;

void coxfit5_b( Sint *maxiter,   Sint *nusedx,    Sint *nvarx, 
	        Sint *strata,    double *beta,    double *u,
	        double *imat2,   double *jmat2,   double *loglik, 
	        Sint *flag,      double *eps,     double *tolerch, 
	        Sint *methodx,   Sint *nfrail,     double *fbeta, 
	        double *fdiag,
                void *fexpr1,    void *fexpr2,     void *rho);

void coxfit5_c(Sint *nusedx,   Sint *nvar,    Sint *strata,
	       Sint *methodx,  double *expect) ;

SEXP coxfit6(SEXP maxiter2,  SEXP time2,   SEXP status2, 
	     SEXP covar2,    SEXP offset2, SEXP weights2,
	     SEXP strata2,   SEXP method2, SEXP eps2, 
	     SEXP toler2,    SEXP ibeta,    SEXP doscale2) ;

void coxmart(Sint   *sn,     Sint   *method,    double *time, 
	     Sint   *status, Sint   * strata,   double *score, 
	     double *wt,     double *expect);

void coxmart2(Sint   *sn,     double *time, 
	     Sint   *status, Sint   * strata,   double *score, 
	     double *wt,     double *resid);

void coxph_wtest(Sint *nvar2, Sint *ntest, double *var, double *b,
                 double *scratch, double *tolerch);

void coxscho(Sint   *nusedx,    Sint   *nvarx,    double *y, 
	     double *covar2,    double *score,    Sint   *strata,  
	     Sint   *method2,   double *work);

void coxscore(Sint   *nx,      Sint   *nvarx,    double *y, 
	      double *covar2,  Sint   *strata,   double *score, 
	      double *weights, Sint   *method,   double *resid2,
	      double *scratch);

double coxsafe(double x);
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

SEXP pyears3b(SEXP   death2,    SEXP   efac2,   SEXP edims2,
	      SEXP   ecut2,     SEXP   expect2, SEXP grpx2,
	      SEXP   x2, 	SEXP   y2,      SEXP times2,
	      SEXP   ngrp2);

double pystep(int nc,        int  *index,  int  *index2,   double *wt, 
	      double *data,  Sint *fac,    Sint *dims,     double **cuts, 
	      double step,   int  edge);

void survdiff2(Sint   *nn,     Sint   *nngroup,    Sint   *nstrat, 
	       double *rho,    double *time,       Sint   *status, 
	       Sint   *group,  Sint   *strata,	   double *obs, 
	       double *exp,    double *var,        double *risk, 
	       double *kaplan);

void survfit4(Sint *n,	Sint *dd,  double *x1,  double *x2) ;

 SEXP survfitci(SEXP ftime2,  SEXP sort12,  SEXP sort22, SEXP ntime2,
                    SEXP status2, SEXP cstate2, SEXP wt2,  SEXP id2,
                    SEXP p2,  SEXP sefit2) ;
  
SEXP survreg6(SEXP maxiter2,   SEXP nvarx,  SEXP y,
	      SEXP ny2,        SEXP covar2, SEXP wtx,
	      SEXP offset2,    SEXP beta2,  SEXP nstratx,
	      SEXP stratax,    SEXP epsx,   SEXP tolx,       
	      SEXP dist,       SEXP expr,   SEXP rho);

SEXP survreg7(SEXP maxiter2,   SEXP nvarx,  SEXP y,
	      SEXP ny2,        SEXP covar2, SEXP wtx,
	      SEXP offset2,    SEXP beta2,  SEXP nstratx,
	      SEXP stratax,    SEXP epsx,   SEXP tolx,
	      SEXP dist,       SEXP dexpr,  SEXP rho,
	      SEXP ptype2,     SEXP pdiag2, SEXP nfrail2,
	      SEXP fgrp2,      SEXP pexpr1, SEXP pexpr2) ;

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

