/*
** Prototypes of all the survival functions
**  Including this in each routine helps prevent mismatched argument errors
*/
SEXP agfit4(SEXP n2,
	    SEXP surv2,      SEXP covar2,    SEXP strata2,
            SEXP weights2,   SEXP offset2,   SEXP ibeta2,
            SEXP sort12,     SEXP sort22,    SEXP method2,
            SEXP maxiter2,   SEXP  eps2,     SEXP tolerance2,
	    SEXP doscale2);


void agfit5a(int  *nusedx,     int  *nvarx,     double *yy, 
	     double *covar2,   double *offset2, double *weights2, 
	     int   *strata,    int    *sort,    double *means,   
	     double *beta,     double *u,       double *loglik, 
	     int  *methodx,    int  *ptype2,    int  *pdiag2,
	     int  *nfrail,     int  *frail2,
	     void *fexpr1,     void *fexpr2,    void *rho,
	     int *docenter) ;

void agfit5b( int  *maxiter,   int  *nusedx,    int  *nvarx, 
	       int  *strata,    double *beta,    double *u,
	       double *imat2,   double *jmat2,   double *loglik, 
	       int  *flag,      double *eps,     double *tolerch, 
	       int  *methodx,   int  *nfrail,     double *fbeta, 
	       double *fdiag,
               void *fexpr1,    void *fexpr2,     void *rho);

void agfit5c( int *nvar);

SEXP agmart3(SEXP nused2,  SEXP surv2,  SEXP score2, SEXP weight2, 
	     SEXP strata2, SEXP sort12, SEXP sort22, SEXP method2);

void agexact(int  *maxiter,  int  *nusedx,   int  *nvarx,   double *start, 
	     double *stop,   int  *event,    double *covar2,double *offset, 
	     int    *strata, double *means,  double *beta,  double *u, 
	     double *imat2,  double loglik[2], int  *flag,  double *work, 
	     int    *work2,  double *eps,    double *tol_chol, double *sctest,
	     int    *nocenter);

void agmart(int    *n,     int    *method,  double *start,   double *stop, 
	    int    *event, double *score,   double *wt,      int    *strata, 
	    double *resid);

SEXP agscore2(SEXP y2,       SEXP covar2,   SEXP strata2, 
	      SEXP score2,   SEXP weights2, SEXP method2);

void agsurv4(int    *ndeath,   double *risk,    double *wt,
             int    *sn,        double *denom,   double *km);

void agsurv5(int  *n2,     int  *nvar2,  int  *dd, double *x1,  
             double *x2,   double *xsum, double *xsum2, 
             double *sum1, double *sum2, double *xbar) ;
SEXP cdecomp(SEXP R2, SEXP time2);

void chinv2  (double **matrix, int n);
int cholesky2(double **matrix, int n, double toler);
void chsolve2(double **matrix, int n, double *y);
void chinv3(double **matrix , int n, int m, double *fdiag);
void chinv5(double **matrix , int n, int flag);
int cholesky3(double **matrix, int n, int m, double *diag, double toler);
int cholesky5(double **matrix, int n, double toler);
void chsolve3(double **matrix, int n, int m, double *diag, double *y);
void chsolve5(double **matrix, int n, double *y, int flag);	

SEXP collapse(SEXP y2,  SEXP x2, SEXP istate2, SEXP id2, SEXP wt2, 
	      SEXP order2) ;
SEXP concordance1(SEXP y, SEXP wt2,  SEXP indx2, SEXP ntree2);
 
SEXP concordance2(SEXP y,     SEXP wt2,  SEXP indx2, SEXP ntree2,
                  SEXP sortstop, SEXP sortstart) ;
SEXP concordance3(SEXP y,        SEXP x2, SEXP wt2, SEXP timewt2, 
                  SEXP sortstop, SEXP doresid2);
SEXP concordance4(SEXP y, SEXP x2, SEXP wt2, SEXP timewt2, 
                  SEXP sortstart, SEXP sortstop, SEXP doresid2); 
SEXP concordance5(SEXP y,        SEXP x2, SEXP wt2, SEXP timewt2, 
                  SEXP sortstop);
SEXP concordance6(SEXP y, SEXP x2, SEXP wt2, SEXP timewt2, 
                  SEXP sortstart, SEXP sortstop); 

void cox_callback(int which, double *coef, double *first, double *second,
		  double *penalty, int *flag, int p, SEXP fexpr, SEXP rho);

SEXP coxcount1(SEXP y2, SEXP strat2) ;
SEXP coxcount2(SEXP y2, SEXP isort1, SEXP isort2, SEXP strat2) ;

void coxdetail(int   *nusedx,    int   *nvarx,     int   *ndeadx, 
	       double *y,        double *covar2,   int   *strata,  
	       double *score,    double *weights,  double *means2, 
	       double *u2,       double *var,      int   *rmat,
	       double *nrisk2,   double *work);
 
SEXP coxexact(SEXP maxiter2,  SEXP y2, 
              SEXP covar2,    SEXP offset2, SEXP strata2,
              SEXP ibeta,     SEXP eps2,    SEXP toler2) ;
 
void coxfit5_a(int *nusedx,     int *nvarx,      double *yy, 
 	       double *covar2,  double *offset2, double *weights2, 
	       int   *strata,   int   *sort,     double *means,   
               double *beta,    double *u,       double *loglik, 
	       int *methodx,    int *ptype2,     int *pdiag2,
	       int *nfrail,     int *frail2,
               void *fexpr1,    void *fexpr2,    void *rho,
	       int  *docenter) ;

void coxfit5_b( int *maxiter,   int *nusedx,     int *nvarx, 
	        int *strata,    double *beta,    double *u,
	        double *imat2,  double *jmat2,   double *loglik, 
	        int *flag,      double *eps,     double *tolerch, 
	        int *methodx,   int *nfrail,     double *fbeta, 
	        double *fdiag,
                void *fexpr1,   void *fexpr2,    void *rho);

void coxfit5_c(int *nusedx,   int *nvar,    int *strata,
	       int *methodx,  double *expect) ;

SEXP coxfit6(SEXP maxiter2,  SEXP time2,   SEXP status2, 
	     SEXP covar2,    SEXP offset2, SEXP weights2,
	     SEXP strata2,   SEXP method2, SEXP eps2, 
	     SEXP toler2,    SEXP ibeta,    SEXP doscale2) ;

void coxmart(int   *sn,     int   *method,    double *time, 
	     int   *status, int   * strata,   double *score, 
	     double *wt,    double *expect);

void coxmart2(int   *sn,     double *time, 
	     int   *status, int   * strata,   double *score, 
	     double *wt,     double *resid);

void coxph_wtest(int *nvar2, int *ntest, double *var, double *b,
                 double *scratch, double *tolerch);

void coxscho(int   *nusedx,    int   *nvarx,    double *y, 
	     double *covar2,    double *score,    int   *strata,  
	     int   *method2,   double *work);

SEXP coxscore2(SEXP y2,       SEXP covar2,   SEXP strata2,
	       SEXP score2,   SEXP weights2, SEXP method2);

double coxsafe(double x);

SEXP coxsurv1(SEXP otime,   SEXP y2,     SEXP weight2, SEXP sort22, 
              SEXP strata2, SEXP xmat2,  SEXP risk2);

SEXP coxsurv2(SEXP otime2, SEXP y2, SEXP weight2,  SEXP sort12, SEXP sort22, 
              SEXP position2,  SEXP strata2, SEXP xmat2, SEXP risk2);

double **dmatrix(double *array, int nrow, int ncol);
int    **imatrix(int *array, int nrow, int ncol);

SEXP finegray(SEXP tstart2, SEXP tstop2,   SEXP ctime2,   SEXP cprob2, 
	      SEXP extend2, SEXP keep2);

SEXP gchol(SEXP matrix2, SEXP toler2);
SEXP gchol_solve(SEXP x2, SEXP y2, SEXP flag2);
SEXP gchol_inv(SEXP matrix, SEXP flag2);

void init_doloop(int min, int max);
int doloop      (int nloops, int *index);

SEXP multicheck(SEXP time12,  SEXP time22, SEXP status2, SEXP id2, 
		SEXP istate2, SEXP sort2);

int *norisk(int n, double *time1, double *time2, double *status, 
	    int *sort1, int *sort2, int *strata);

void pyears1(int   *sn,      int   *sny,      int   *sdoevent, 
	     double *sy,      double *wt,       
	     int   *sedim,   int   *efac, 
	     int   *edims,   double *secut,    double *expect, 
	     double *sedata,  int   *sodim,    int   *ofac, 
	     int   *odims,   double *socut,    int   *smethod, 
	     double *sodata,  double *pyears,   double *pn, 
	     double *pcount,  double *pexpect,  double *offtable);

void pyears2(int   *sn,      int   *sny,   int   *sdoevent, 
	     double *sy,      double *wt,    
	     int   *sodim,   int   *ofac, 
	     int   *odims,   double *socut, double *sodata,
	     double *pyears,  double *pn,    double *pcount, 
	     double *offtable);

SEXP pyears3b(SEXP   death2,    SEXP   efac2,   SEXP edims2,
	      SEXP   ecut2,     SEXP   expect2, SEXP grpx2,
	      SEXP   x2, 	SEXP   y2,      SEXP times2,
	      SEXP   ngrp2);

double pystep(int nc,        int  *index,  int  *index2,   double *wt, 
	      double *data,  int *fac,    int *dims,     double **cuts, 
	      double step,   int  edge);

void survdiff2(int   *nn,     int   *nngroup,    int   *nstrat, 
	       double *rho,    double *time,       int   *status, 
	       int   *group,  int   *strata,	   double *obs, 
	       double *exp,    double *var,        double *risk, 
	       double *kaplan);

void survfit4(int *n,	int *dd,  double *x1,  double *x2) ;

SEXP survfitci(SEXP ftime2,       SEXP sort12,  SEXP sort22, SEXP ntime2,
                    SEXP status2, SEXP cstate2, SEXP wt2,    SEXP id2,
                    SEXP p2,      SEXP i02,     SEXP sefit2) ;
  
SEXP survfitkm(SEXP y2,     SEXP weight2,  SEXP sort12, SEXP sort22, 
               SEXP type2,  SEXP id2,      SEXP nid2,   SEXP position2, 
               SEXP influence2) ;

SEXP survfitresid(SEXP Y2,      SEXP sort12,  SEXP sort22,  SEXP cstate2, 
		  SEXP wt2,     SEXP p02,     SEXP i02,     SEXP otime2,  
		  SEXP starttime2, SEXP doauc2);
 
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
		 double *beta,   int dist,     int *strat,     double *offset,
		 double *time1,  double *time2, double *status, double *wt,
		 double **covar, double **imat, double **JJ,    double *u, 
		 SEXP expr,      SEXP rho,      double *dummy,  int nf,
		 int *frail,    double *fdiag, double *jdiag );

double survregc2(int n,          int nvar,     int nstrat,      int whichcase,
		 double *beta,   int dist,     int *strat,     double *offset,
		 double *time1,  double *time2, double *status, double *wt,
		 double **covar, double **imat, double **JJ,    double *u, 
		 SEXP expr,      SEXP rho,      double *dummy,  int nf,
		 int *frail,    double *fdiag, double *jdiag );

void survpenal(int whichcase, int nfrail,    int  nvar2,    double **hmat, 
	       double **JJ,   double *hdiag, double *jdiag,
	       double *u,     double *beta,  double *loglik,
	       int ptype,     int pdiag,     SEXP pexpr1,   double *cptr1, 
	       SEXP pexpr2,   double *cptr2, SEXP rho);

SEXP survsplit(SEXP tstart2,  SEXP tstop2,  SEXP cut2);

SEXP tmerge (SEXP id2,  SEXP time1x, SEXP newx2,
	     SEXP nid2, SEXP ntime2, SEXP x2,  SEXP indx2); 
SEXP tmerge2(SEXP id2,  SEXP time1x, SEXP nid2, SEXP ntime2);
SEXP tmerge3(SEXP id2, SEXP miss2);

void walkup(double *nwt, double* twt, int index, double sums[3], int ntree);
void addin(double *nwt, double *twt, int index, double wt);

SEXP zph1(SEXP gt2,    SEXP y2, 
	  SEXP covar2, SEXP eta2,  SEXP weights2,
	  SEXP strata2,SEXP method2, SEXP sort2);

SEXP zph2(SEXP gt2,    SEXP y2, 
	  SEXP covar2, SEXP eta2,  SEXP weights2,
	  SEXP strata2,SEXP method2, SEXP sort12, SEXP sort22);
