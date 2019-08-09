/*
** This file causes the entry points of my .C routines to be preloaded
** Added at the request of R-core.
** It adds one more layer of protection by declaring the number of arguments,
**  and perhaps a tiny bit of speed
*/
#include "survS.h"
#include "R_ext/Rdynload.h"
#include "Rversion.h"
#include "survproto.h"

static const R_CMethodDef Centries[] = {
    {"Cagfit5a",    (DL_FUNC) &agfit5a,  20},
    {"Cagfit5b",    (DL_FUNC) &agfit5b,  19},
    {"Cagfit5c",    (DL_FUNC) &agfit5c,   1},
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
    {"Cagfit4",       (DL_FUNC) &agfit4,      13},
    {"Cagmart3",      (DL_FUNC) &agmart3,      6},
    {"Cconcordance1", (DL_FUNC) &concordance1, 4}, 
    {"Cconcordance2", (DL_FUNC) &concordance2, 6}, 
    {"Cconcordance3", (DL_FUNC) &concordance3, 6},
    {"Cconcordance4", (DL_FUNC) &concordance4, 7},
    {"Ccoxcount1",    (DL_FUNC) &coxcount1,    2},
    {"Ccoxcount2",    (DL_FUNC) &coxcount2,    4},
    {"Ccoxexact",     (DL_FUNC) &coxexact,     8},
    {"Ccoxfit6",      (DL_FUNC) &coxfit6,     12},
    {"Ccoxsurv1",     (DL_FUNC) &coxsurv1,    11},
    {"Ccoxsurv2",     (DL_FUNC) &coxsurv2,     9},
    {"Cfinegray",     (DL_FUNC) &finegray,     6},
    {"Cgchol",        (DL_FUNC) &gchol,        2},
    {"Cgchol_solve",  (DL_FUNC) &gchol_solve,  3},
    {"Cgchol_inv",    (DL_FUNC) &gchol_inv,    2},
    {"Cmulticheck",   (DL_FUNC) &multicheck,   5},
    {"Cpyears3b",     (DL_FUNC) &pyears3b,    10},
    {"Csurvfitci",    (DL_FUNC) &survfitci,   11},
    {"Csurvfitkm",    (DL_FUNC) &survfitkm,    9},
    {"Csurvreg6",     (DL_FUNC) &survreg6,    15},
    {"Csurvreg7",     (DL_FUNC) &survreg7,    21},
    {"Csurvsplit",    (DL_FUNC) &survsplit,    3},
    {"Ctmerge",       (DL_FUNC) &tmerge,       7},
    {"Czph1",         (DL_FUNC) &zph1,         8},
    {"Czph2",         (DL_FUNC) &zph2,         9},
    {NULL, NULL, 0}
};

void R_init_survival(DllInfo *dll){
    R_registerRoutines(dll, Centries, Callentries, NULL, NULL);

    /* The following line makes only those routines defined above
       available to outside packages, i.e., internal things like
       dmatrix() are now invisible.
    */
    R_useDynamicSymbols(dll, FALSE); 
    /*
    ** This line makes them only be available via the symbols above
    **  i.e., .Call("tmerge", ) won't work but .Call(Ctmerge, )  will
    ** This feature was added in version 2.16
    */
#if defined(R_VERSION) && R_VERSION >= R_Version(2, 16, 0)
    R_forceSymbols(dll, TRUE);
#endif
}
    
