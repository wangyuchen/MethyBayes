#include <R_ext/Rdynload.h>
#include <R.h>

void F77_NAME(mcmc)(int *l, int *g0, int *g1, int *g, double *p, int *result);

static const R_FortranMethodDef fMethods[] = {
    {"mcmc_froutine",  (DL_FUNC) &mcmc_, 6},
    {NULL, NULL, 0}
};

void R_init_MethyBayes(DllInfo *info)
{
    R_registerRoutines(info, NULL, NULL, fMethods, NULL);
}
