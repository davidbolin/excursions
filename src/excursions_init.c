#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <Rinternals.h>

/* .C calls */
extern void reordering(void *, void *, void *, void *, void *);
extern void shapeInt(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void testRand(void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"reordering", (DL_FUNC) &reordering,  5},
    {"shapeInt",   (DL_FUNC) &shapeInt,   10},
    {"testRand",   (DL_FUNC) &testRand,    3},
    {NULL, NULL, 0}
};

/* .Call calls */
extern SEXP Qinv(SEXP, SEXP, SEXP, SEXP);
extern SEXP excursions_openmp_info(void);

static const R_CallMethodDef CallEntries[] = {
    {"Qinv", (DL_FUNC) &Qinv, 4},
    {"excursions_openmp_info", (DL_FUNC) &excursions_openmp_info, 0},
    {NULL, NULL, 0}
};

void R_init_excursions(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
