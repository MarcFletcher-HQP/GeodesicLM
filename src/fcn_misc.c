#ifdef LOCAL_HEADERS
#include "R.h"
#include "Rdefines.h"
#include "GeodesicLM.h"

#else
#include <R.h>
#include <Rdefines.h>
#include "GeodesicLM.h"

#endif

void fcn_ja(int *ldfjac, int *n, double *par, double **fjac) { Rprintf("Jacobians!\n"); }
void fcn_ac(int *ldfacc, int *n, double *par, double *v, double *facc) { Rprintf("Acceleration!\n"); }