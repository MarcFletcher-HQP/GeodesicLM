#ifdef LOCAL_HEADERS
#include "R.h"
#include "Rdefines.h"
#include "GeodesicLM.h"

#else
#include <R.h>
#include <Rdefines.h>
#include "GeodesicLM.h"

#endif

/* This function is passed to lmdif.f, the function is evaluated at 'par' and the result 
is stored in 'fvec'. */
void fcn_lm(int *m, int *n, double *par, double *fvec)
{
    int i;
	SEXP sexp_fvec;

    /* Rprintf("fcn-lmdif calling...\n"); */
	// Check that 'par' is within the bounds set by 'lower' and 'upper'
    if (IS_NUMERIC(OS->par))
      for (i = 0; i < *n; i++) {
		if(par[i] < NUMERIC_POINTER(OS->lower)[i])
			par[i] = NUMERIC_POINTER(OS->lower)[i];
		if(par[i] > NUMERIC_POINTER(OS->upper)[i])
			par[i] = NUMERIC_POINTER(OS->upper)[i];
		Rprintf("%f ", par[i]);
		NUMERIC_POINTER(OS->par)[i] = par[i];
      }
    else
      for (i = 0; i < *n; i++) {
		if(par[i] < NUMERIC_POINTER(OS->lower)[i])
			par[i] = NUMERIC_POINTER(OS->lower)[i];
		if(par[i] > NUMERIC_POINTER(OS->upper)[i])
			par[i] = NUMERIC_POINTER(OS->upper)[i];
		Rprintf("%f", par[i]);
		NUMERIC_POINTER(VECTOR_ELT(OS->par, i))[0] = par[i];
      }
	Rprintf("\n");

	Rprintf("Updating fvec!\n");
	SETCADR(OS->fcall, OS->par);
	PROTECT(sexp_fvec = eval(OS->fcall, OS->env));

	Rprintf("fvec = ");
	for (i = 0; i < *m; i++) {
		fvec[i] = NUMERIC_POINTER(sexp_fvec)[i];
		Rprintf("%g ", fvec[i]);
	}
	Rprintf("\n");
	UNPROTECT(1);
}
