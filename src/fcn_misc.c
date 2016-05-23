#ifdef LOCAL_HEADERS
#include "R.h"
#include "Rdefines.h"
#include "GeodesicLM.h"

#else
#include <R.h>
#include <Rdefines.h>
#include "GeodesicLM.h"

#endif

void fcn_ja(int *ldfjac, int *n, double *par, double *fjac) {
	
	int i;
	SEXP sexp_fjac;

	/* Rprintf("fcn-lmdif calling...\n"); */
	// Update value of 'par' stored in OS

	SETCADR(OS->jcall, OS->par);
	PROTECT(sexp_fjac = eval(OS->jcall, OS->env));

	for (i = 0; i < *m; i++)
		fjac[i] = NUMERIC_POINTER(sexp_jac)[i];
	UNPROTECT(1);
}
void fcn_ac(int *ldfacc, int *n, double *par, double *v, double *facc) {

	int i;
	SEXP sexp_facc;

	/* Rprintf("fcn-lmdif calling...\n"); */
	// Update value of 'par' stored in OS

	SETCADR(OS->acall, OS->par);
	PROTECT(sexp_facc = eval(OS->acall, OS->env));

	for (i = 0; i < *m; i++)
		facc[i] = NUMERIC_POINTER(sexp_acc)[i];
	UNPROTECT(1);
}