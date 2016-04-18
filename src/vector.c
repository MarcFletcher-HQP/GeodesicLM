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

int *int_vector(int n)
{
    return (int*) R_alloc(n, sizeof(int));
}
