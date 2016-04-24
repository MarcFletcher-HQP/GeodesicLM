#ifdef LOCAL_HEADERS
#include "R.h"
#include "Rdefines.h"
#include "GeodesicLM.h"

#else
#include <R.h>
#include <Rdefines.h>
#include "GeodesicLM.h"

#endif

/* This function is called by geodesiclm.f90 after each iteration, the function is evaluated at 'par' and the result
is stored in 'fvec'. */
void fcn_call(int *m, int *n, double *par, double *v, double *a, double *fvec, double *fjac,
	double *acc, double *lam, double *dtd, double *fvec_new, int *accepted, int *info)
{
	// Set iflag if niter reaches the maximum
	if (OS->niter == OS->maxiter)
		*info = -1;
}