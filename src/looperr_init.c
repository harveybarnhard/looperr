#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _looperr_fastols(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _looperr_fastols_by(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _looperr_fastols_bywt(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _looperr_loclin_diffX(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _looperr_loclin_diffX_unif(SEXP, SEXP, SEXP, SEXP);
extern SEXP _looperr_loclin_sameX(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _looperr_loclin_sameX_unif(SEXP, SEXP, SEXP);
extern SEXP _looperr_loclin_sameX_unif_by(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_looperr_fastols",              (DL_FUNC) &_looperr_fastols,              5},
    {"_looperr_fastols_by",           (DL_FUNC) &_looperr_fastols_by,           6},
    {"_looperr_fastols_bywt",         (DL_FUNC) &_looperr_fastols_bywt,         7},
    {"_looperr_loclin_diffX",         (DL_FUNC) &_looperr_loclin_diffX,         6},
    {"_looperr_loclin_diffX_unif",    (DL_FUNC) &_looperr_loclin_diffX_unif,    4},
    {"_looperr_loclin_sameX",         (DL_FUNC) &_looperr_loclin_sameX,         5},
    {"_looperr_loclin_sameX_unif",    (DL_FUNC) &_looperr_loclin_sameX_unif,    3},
    {"_looperr_loclin_sameX_unif_by", (DL_FUNC) &_looperr_loclin_sameX_unif_by, 5},
    {NULL, NULL, 0}
};

void R_init_looperr(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
