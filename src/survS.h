/* 
**  This file started out with support for Splus, then morphed to allow
**  either R or Splus (based on ifdef lines), and now is R only.
*/
#include "R.h"
#include "Rinternals.h"
#include <R_ext/Utils.h>  

/*
** Memory defined with ALLOC is removed automatically by S.
**  That with "Calloc" I have to remove myself.  Use the
**  latter for objects that need to to persist between calls.
*/
#define ALLOC(a,b)  R_alloc(a,b)
#define CALLOC(a,b) R_Calloc(a,b)
#define FREE(a)     R_Free(a)
