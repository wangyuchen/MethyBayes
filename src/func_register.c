#include "MethyBayes.h"
#include <R_ext/Rdynload.h>


void C_mcmc(int *l, int *g0, int *g1, int *g, double *p, int *result) {
    F77_CALL(mcmc)(l, g0, g1, g, p, result);
}


static const R_CMethodDef cMethods[] = {
    {"C_mcmc",  (DL_FUNC) &C_mcmc, 6},
    NULL
};

void R_init_MethyBayes(DllInfo *info)
{
    /* Register the .C and .Call routines.
    No .Fortran() or .External() routines,
    so pass those arrays as NULL.
    */
    R_registerRoutines(info, cMethods, NULL, NULL, NULL);
}
