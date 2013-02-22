/*
** This file causes the entry points of my .C routines to be preloaded
** Added at the request of R-core.
** It adds one more layer of protection by declaring the number of arguments,
**  and perhaps a tiny bit of speed
*/
#include "survS.h"
#include "R_ext/Rdynload.h"
#include "survproto.h"

static const R_CMethodDef Centries[] = {
    {"Cagfit3",     (DL_FUNC) &agfit3,    23},
    {"Cagfit5a",    (DL_FUNC) &agfit5_a,  20},
    {"Cagfit5b",    (DL_FUNC) &agfit5_b,  19},
    {"Cagfit5c",    (DL_FUNC) &agfit5_c,   5},
    {"Cagsurv3",    (DL_FUNC) &agsurv3,   19},
    {"Cagsurv4",    (DL_FUNC) &agsurv4,    6},
    {"Cagsurv5",    (DL_FUNC) &agsurv5,   10},
    {"Cagexact",    (DL_FUNC) &agexact,   20},
    {"Cagmart",     (DL_FUNC) &agmart,     9},
    {"Cagmart2",    (DL_FUNC) &agmart2,   13},
    {"Cagscore",    (DL_FUNC) &agscore,   10},
    {"Ccoxdetail",  (DL_FUNC) &coxdetail, 14},
    {"Ccoxfit5a",   (DL_FUNC) &coxfit5_a, 20},
    {"Ccoxfit5b",   (DL_FUNC) &coxfit5_b, 19},
    {"Ccoxfit5c",   (DL_FUNC) &coxfit5_c,  5},
    {"Ccoxmart",    (DL_FUNC) &coxmart,    8},
    {"Ccoxmart2",   (DL_FUNC) &coxmart2,   7},
    {"Ccoxph_wtest",(DL_FUNC) &coxph_wtest,6},
    {"Ccoxscho",    (DL_FUNC) &coxscho,    8},
    {"Ccoxscore",   (DL_FUNC) &coxscore,  10},
    {"Cpyears1",    (DL_FUNC) &pyears1,   22},
    {"Cpyears2",    (DL_FUNC) &pyears2,   14},
    {"Csurvdiff2",  (DL_FUNC) &survdiff2, 13},
    {"Csurvfit4",   (DL_FUNC) &survfit4,   4},
    {NULL, NULL, 0}
};

static const R_CallMethodDef Callentries[] = {
    {"Cconcordance1", (DL_FUNC) &concordance1, 4}, 
    {"Cconcordance2", (DL_FUNC) &concordance2, 6}, 
    {"Ccoxcount1",    (DL_FUNC) &coxcount1,    2},
    {"Ccoxcount2",    (DL_FUNC) &coxcount2,    4},
    {"Ccoxexact",     (DL_FUNC) &coxexact,     8},
    {"Ccoxfit6",      (DL_FUNC) &coxfit6,     12},
    {"Cpyears3b",     (DL_FUNC) &pyears3b,    10},
    {"Csurvfitci",    (DL_FUNC) &survfitci,   10},
    {"Csurvreg6",     (DL_FUNC) &survreg6,    15},
    {"Csurvreg7",     (DL_FUNC) &survreg7,    21},
    {NULL, NULL, 0}
};

void R_init_survival(DllInfo *dll){
    R_registerRoutines(dll, Centries, Callentries, NULL, NULL);

    /* My take on the documentation is that adding the following line
       will make symbols available ONLY through the above tables.
       Anyone who then tried to link to my C code would be SOL.
       It also wouldn't work with .C(routines[1], ....
    */
   /* R_useDynamicSymbols(dll, FALSE);  */
}
    
