#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void Qinv(void *, void *, void *, void *, void *, void *);
extern void reordering(void *, void *, void *, void *, void *);
extern void shapeInt(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void testRand(void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"Qinv",       (DL_FUNC) &Qinv,        6},
    {"reordering", (DL_FUNC) &reordering,  5},
    {"shapeInt",   (DL_FUNC) &shapeInt,   10},
    {"testRand",   (DL_FUNC) &testRand,    3},
    {NULL, NULL, 0}
};

void R_init_excursions(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
