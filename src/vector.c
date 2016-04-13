#ifdef LOCAL_HEADERS
#include "R.h"
#include "Rdefines.h"

#else
#include <R.h>
#include <Rdefines.h>

#endif

double *real_vector(int n)
{
    return (double*) R_alloc(n, sizeof(double));
}

// See the project: 'matrix_tests' for a case where this works
double **real_matrix(int m, int n)
{
	int		i, j;
	double	**A;

	A = malloc(m * sizeof(double*));
	for (i = 0; i < m; i++) {
		A[i] = malloc(n * sizeof(double*));
		for (j = 0; j < n; j++) {
			A[i][j] = 0;
		}
	}

	return A;
}

int *int_vector(int n)
{
    return (int*) R_alloc(n, sizeof(int));
}

// Need to test this in matrix_tests
SEXP SET_ARRAY_DIM(SEXP x, int m, int n) 
{
	SEXP	out_dim, out;
	
	PROTECT(out = duplicate(x));
	PROTECT(out_dim = allocVector(INTSXP, 2));

	INTEGER(out_dim)[0] = m;
	INTEGER(out_dim)[1] = n;

	setAttrib(out, R_DimSymbol, out_dim);

	UNPROTECT(2);
	return out;
}
