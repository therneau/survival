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
    {"Cagfit5a",    (DL_FUNC) &agfit5a,  21},
    {"Cagfit5b",    (DL_FUNC) &agfit5b,  19},
    {"Cagfit5c",    (DL_FUNC) &agfit5c,   1},
    {"Cagsurv4",    (DL_FUNC) &agsurv4,    6},
    {"Cagsurv5",    (DL_FUNC) &agsurv5,   10},
    {"Cagexact",    (DL_FUNC) &agexact,   21},
    {"Cagmart",     (DL_FUNC) &agmart,     9},
    {"Ccoxdetail",  (DL_FUNC) &coxdetail, 15},
    {"Ccoxfit5a",   (DL_FUNC) &coxfit5_a, 21},
    {"Ccoxfit5b",   (DL_FUNC) &coxfit5_b, 19},
    {"Ccoxfit5c",   (DL_FUNC) &coxfit5_c,  5},
    {"Ccoxmart",    (DL_FUNC) &coxmart,    8},
    {"Ccoxmart2",   (DL_FUNC) &coxmart2,   7},
    {"Ccoxph_wtest",(DL_FUNC) &coxph_wtest,6},
    {"Ccoxscho",    (DL_FUNC) &coxscho,    8},
    {"Cpyears1",    (DL_FUNC) &pyears1,   22},
    {"Cpyears2",    (DL_FUNC) &pyears2,   14},
    {"Csurvdiff2",  (DL_FUNC) &survdiff2, 13},
    {"Csurvfit4",   (DL_FUNC) &survfit4,   4},
    {NULL, NULL, 0}
};

static const R_CallMethodDef Callentries[] = {
    {"Cagfit4",       (DL_FUNC) &agfit4,      14},
    {"Cagmart3",      (DL_FUNC) &agmart3,      8},
    {"Cagscore2",     (DL_FUNC) &agscore2,     6},
    {"Cagscore3",     (DL_FUNC) &agscore3,     7},
    {"Ccdecomp",      (DL_FUNC) &cdecomp,      2},
    {"Ccollapse",     (DL_FUNC) &collapse,     6},
    {"Cconcordance1", (DL_FUNC) &concordance1, 4}, 
    {"Cconcordance2", (DL_FUNC) &concordance2, 6}, 
    {"Cconcordance3", (DL_FUNC) &concordance3, 6},
    {"Cconcordance4", (DL_FUNC) &concordance4, 7},
    {"Cconcordance5", (DL_FUNC) &concordance5, 5},
    {"Cconcordance6", (DL_FUNC) &concordance6, 6},
    {"Ccoxcount1",    (DL_FUNC) &coxcount1,    2},
    {"Ccoxcount2",    (DL_FUNC) &coxcount2,    4},
    {"Ccoxexact",     (DL_FUNC) &coxexact,     8},
    {"Ccoxfit6",      (DL_FUNC) &coxfit6,     12},
    {"Ccoxscore2",    (DL_FUNC) &coxscore2,    6},
    {"Ccoxsurv1",     (DL_FUNC) &coxsurv1,     7},
    {"Ccoxsurv2",     (DL_FUNC) &coxsurv2,     9},
    {"Ccoxsurv3",     (DL_FUNC) &coxsurv3,     7},
    {"Ccoxsurv4",     (DL_FUNC) &coxsurv4,     8},
    {"Cfastkm1",      (DL_FUNC) &fastkm1,      3},
    {"Cfastkm2",      (DL_FUNC) &fastkm1,      4},
    {"Cfinegray",     (DL_FUNC) &finegray,     6},
    {"Cgchol",        (DL_FUNC) &gchol,        2},
    {"Cgchol_solve",  (DL_FUNC) &gchol_solve,  3},
    {"Cgchol_inv",    (DL_FUNC) &gchol_inv,    2},
    {"Cmulticheck",   (DL_FUNC) &multicheck,   6},
    {"Cpyears3b",     (DL_FUNC) &pyears3b,    10},
    {"Cresidcsum",    (DL_FUNC) &residcsum,    2},
    {"Csurvfitaj",    (DL_FUNC) &survfitaj,   16},
    {"Csurvfitkm",    (DL_FUNC) &survfitkm,   11},
    {"Csurvfitresid", (DL_FUNC) &survfitresid,10},
    {"Csurvreg6",     (DL_FUNC) &survreg6,    15},
    {"Csurvreg7",     (DL_FUNC) &survreg7,    21},
    {"Csurvsplit",    (DL_FUNC) &survsplit,    3},
    {"Ctmerge",       (DL_FUNC) &tmerge,       7},
    {"Ctmerge2",      (DL_FUNC) &tmerge2,      4},
    {"Ctmerge3",      (DL_FUNC) &tmerge3,      2},
    {"Ctwoclust",     (DL_FUNC) &twoclust,     3},
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
    **  i.e., .Call("tmerge", ) won't work but .Call(Ctmerge, ) will.
    ** This feature was added in version 2.16
    */
#if defined(R_VERSION) && R_VERSION >= R_Version(2, 16, 0)
    R_forceSymbols(dll, TRUE);
#endif
}
    
