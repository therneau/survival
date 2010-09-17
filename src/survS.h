/* SCCS $Id: survS.h 11252 2009-03-19 13:46:26Z tlumley $ 
**
**   The S.h file defines a few things that I need, and hundreds that I don't.
** In particular, on some architectures, it defines a variable "time"
** which of course conflicts with lots of my C-code, 'time' being a natural
** variable name for survival models.
**   Thanks to Brian Ripley for suggesting a machine independent way of
** fixing this.
**
** The S_alloc function changed it's argument list from Splus version
**   4 to 5, and there is a different one for R.
**   The ALLOC macro allows me to have common C code for all versions,
**   with only this file "survS.h" changed.
*/


#include "R.h"
#include "Rinternals.h"
#ifdef USING_R
   /* typedef int Sint; */
#define S_EVALUATOR    /* Turn this into a "blank line" in R */
#else
/*
** Splus definitions, to use R type calls
*/
typedef long Sint;
/* 
**  At this point in time (Splus 8.0.1) I need to add a little
**   to the Insightful definitions.  (They are in the process
**   of improving Rinternals, so this may well go away.)  The 
**   two functions below are defined as "not supported".  I need
**   only certain cases of defineVar and eval, so can safely map them.
**   I am using the 8.1 R*.h files courtesy of Bill Dunlap
*/
#ifdef  defineVar
#undef  defineVar
#endif
#define defineVar(a,b,c) ASSIGN_IN_FRAME(a,b, INTEGER_VALUE(c))

#ifdef  eval
#undef  eval
#endif
#define eval(a, b)  EVAL_IN_FRAME(a, INTEGER_VALUE(b))

/*
** These two refer to undefined functions, so use the 8.0.1 defs
*/
#ifdef asInteger
#undef asInteger
#endif
#define asInteger(a) INTEGER_VALUE(a)
#ifdef asReal
#undef asReal
#endif
#define asReal(a) NUMERIC_VALUE(a)

#endif

/*
** Memory defined with ALLOC is removed automatically by S.
**  That with "Calloc" I have to remove myself.  Use the
**  latter for objects that need to to persist between calls.
*/
#ifdef USING_R
#define ALLOC(a,b)  R_alloc(a,b)
#else
#define ALLOC(a,b)  S_alloc(a,b)
#endif

/*
** Prototype for callback function
**
*/
#ifdef USING_R
void cox_callback(int which, double *coef, double *first, double *second,
		  double *penalty, int *flag, int p, SEXP fexpr, SEXP rho);
#endif
