/*
** callback routines for the coxph frailty interface
*** This is completely  rewritten for R (TSL April 99)
*/
#include "survS.h"
#include "Rinternals.h"
#include "localization.h"
/*
** This part is called by the coxfit4 function, to get the penalty terms
*/

/**
*** We need to call fexpr in environment rho in R
*** These return lists, which need to be processed to get at the elements
***
*** Also, we pass in the length of coef rather than messing about in
*** the calling environment to find it.
**/

void cox_callback (int which,      double *coef,    double *first, 
		   double *second, double *penalty, int *flag,
		   int p,          SEXP fexpr,      SEXP rho) {
    SEXP coxlist, temp,data,index;
    int i;

    /** copy coef into R vector */

    PROTECT(data=allocVector(REALSXP,p));
    for (i=0;i<p;i++){
      REAL(data)[i]=coef[i];
	}

    /** eval function */
    PROTECT(temp=lang2(fexpr, data));
    PROTECT(coxlist=eval(temp,rho));
    UNPROTECT(3);
    PROTECT(coxlist);
    /* stick it back in the calling frame */
    if (which==1)
      setVar(install("coxlist1"),coxlist,rho);
    else
      setVar(install("coxlist2"),coxlist,rho);
   /* Grab the updated values from the list */
    PROTECT(index=mkString("coef"));
    PROTECT(temp=lang3(install("[["),coxlist,index));
    PROTECT(data=eval(temp,rho));
    if (!isNumeric(data))
                error(_("coef: invalid type\n"));
    for (i=0;i<length(data);i++){
      coef[i]=REAL(data)[i];
    }
    UNPROTECT(3);
    PROTECT(index=mkString("first"));
    PROTECT(temp=lang3(install("[["),coxlist,index));
    PROTECT(data=eval(temp,rho));
    if (!isNumeric(data))
                error(_("first: invalid type\n"));
    for (i=0;i<length(data);i++){
      first[i]=REAL(data)[i];
      /* printf("%g,",first[i]);*/
    }
    UNPROTECT(3);
    PROTECT(index=mkString("second"));
    PROTECT(temp=lang3(install("[["),coxlist,index));
    PROTECT(data=eval(temp,rho));
    if (!isNumeric(data))
                error(_("second: invalid type\n"));
    for (i=0;i<length(data);i++){
      second[i]=REAL(data)[i];
    }
    UNPROTECT(3);
    PROTECT(index=mkString("flag"));
    PROTECT(temp=lang3(install("[["),coxlist,index));
    PROTECT(data=eval(temp,rho));
    if (!(isInteger(data) | isLogical(data)))
                error(_("flag: invalid type\n"));
    for (i=0;i<length(data);i++){
      flag[i]=LOGICAL(data)[i];
    }
    UNPROTECT(3);
    PROTECT(index=mkString("penalty"));
    PROTECT(temp=lang3(install("[["),coxlist,index));
    PROTECT(data=eval(temp,rho));
    if (!isNumeric(data))
                error(_("penalty: invalid type\n"));
    for (i=0;i<length(data);i++){
      penalty[i]=REAL(data)[i];
    }
    UNPROTECT(3);
    /* clean up */
    UNPROTECT(1); /*coxlist*/

}










