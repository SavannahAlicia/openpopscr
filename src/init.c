#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _openpopscr_C_calc_D(SEXP, SEXP, SEXP, SEXP);
extern SEXP _openpopscr_C_calc_llk(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _openpopscr_C_calc_move_llk(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _openpopscr_C_calc_pdet(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _openpopscr_C_calc_move_pdet(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _openpopscr_C_calc_pr_capture(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _openpopscr_C_calc_alpha(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _openpopscr_C_calc_beta(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _openpopscr_C_calc_movealpha(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _openpopscr_C_calc_movebeta(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_openpopscr_C_calc_D",          (DL_FUNC) &_openpopscr_C_calc_D,          4},
    {"_openpopscr_C_calc_llk",        (DL_FUNC) &_openpopscr_C_calc_llk,        8},
    {"_openpopscr_C_calc_move_llk",   (DL_FUNC) &_openpopscr_C_calc_move_llk,   16},
    {"_openpopscr_C_calc_pdet",       (DL_FUNC) &_openpopscr_C_calc_pdet,       5},
    {"_openpopscr_C_calc_move_pdet",  (DL_FUNC) &_openpopscr_C_calc_move_pdet,  14},
    {"_openpopscr_C_calc_pr_capture", (DL_FUNC) &_openpopscr_C_calc_pr_capture, 17},
    {"_openpopscr_C_calc_alpha",      (DL_FUNC) &_openpopscr_C_calc_alpha,      8},
    {"_openpopscr_C_calc_beta",      (DL_FUNC) &_openpopscr_C_calc_beta,        8},
    {"_openpopscr_C_calc_movealpha",      (DL_FUNC) &_openpopscr_C_calc_movealpha,  16},
    {"_openpopscr_C_calc_movebeta",      (DL_FUNC) &_openpopscr_C_calc_movebeta,    16},
    {NULL, NULL, 0}
};

void R_init_openpopscr(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
