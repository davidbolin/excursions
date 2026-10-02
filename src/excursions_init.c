#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <Rinternals.h>

/* .C calls */
extern void reordering(void *, void *, void *, void *, void *);
extern void shapeInt(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void shapeIntProbe(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void testRand(void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"reordering", (DL_FUNC) &reordering,  5},
    {"shapeInt",   (DL_FUNC) &shapeInt,   10},
    {"shapeIntProbe", (DL_FUNC) &shapeIntProbe, 17},
    {"testRand",   (DL_FUNC) &testRand,    3},
    {NULL, NULL, 0}
};

/* .Call calls */
extern SEXP Qinv(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP Qinv_selected(SEXP, SEXP, SEXP, SEXP);
extern SEXP excursions_openmp_info(void);
extern SEXP shapeIntCall(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP regions_bvn_lower(SEXP, SEXP, SEXP);
extern SEXP regions_components(SEXP, SEXP, SEXP);
extern SEXP regions_grow(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP regions_prominence(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"Qinv", (DL_FUNC) &Qinv, 5},
    {"Qinv_selected", (DL_FUNC) &Qinv_selected, 4},
    {"excursions_openmp_info", (DL_FUNC) &excursions_openmp_info, 0},
    {"shapeIntCall", (DL_FUNC) &shapeIntCall, 9},
    {"regions_bvn_lower", (DL_FUNC) &regions_bvn_lower, 3},
    {"regions_components", (DL_FUNC) &regions_components, 3},
    {"regions_grow", (DL_FUNC) &regions_grow, 8},
    {"regions_prominence", (DL_FUNC) &regions_prominence, 4},
    {NULL, NULL, 0}
};

void R_init_excursions(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
