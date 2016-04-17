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
void fcn_call(int *m, int *n, double *par, double *v, double *a, double *fvec, double **fjac,
	double *acc, double *lam, double **dtd, double *fvec_new, int *accepted, int *info)
{
	int i;
	double sumf;
	SEXP sexp_fvec;

	/* Note, previously this section would only run if a variable called 'iflag' was set to zero. There is no flag passed as an argument in geodesicLM,
	as such the flag has been moved into opt_struct and must be initialised to zero before the call to geodesicLM */
	if (OS->converged == 0) {
		if (OS->nprint > 0) {
			Rprintf("It. %4d, RSS = %10g, Par. =", OS->niter,
				OS->rsstrace[OS->niter]);
			for (i = 0; i < *n; i++)
				Rprintf(" % 10g", par[i]);
			Rprintf("\n");
		}
		OS->niter++;	// increment the number of iterations
	}
	// If algorithm has finished then finalise fvec values and trace.
	else if (OS->converged == 1 || OS->converged == 2 || OS->converged == 3 || OS->converged == 4 ||
		OS->converged == 5 || OS->converged == 6 || OS->converged == 7) {
		SETCADR(OS->fcall, OS->par);
		PROTECT(sexp_fvec = eval(OS->fcall, OS->env));
		Rprintf("OS->converged = %i", OS->converged);
		for (i = 0; i < *m; i++)
			fvec[i] = NUMERIC_POINTER(sexp_fvec)[i];
		UNPROTECT(1);
		sumf = 0;
		for (i = 0; i < *m; i++)
			sumf += (fvec[i] * fvec[i]);

		OS->rsstrace[OS->niter] = sumf;
	}

	if (OS->niter == OS->maxiter)
		*info = -1;
}