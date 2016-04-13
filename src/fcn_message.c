#ifdef LOCAL_HEADERS
#include "R.h"
#include "Rdefines.h"
#include "GeodesicLM.h"

#else
#include <R.h>
#include <Rdefines.h>
#include "GeodesicLM.h"

#endif

char *fcn_message(char *msg, int converged, int info, int n, int nit, int nfev, int njev, int naev)
{
	Rprintf("converged = %i \n info = %i \n", converged, info);
    if      (converged == 1)
        sprintf(msg, "cosine of the angle between the residual vector and the range of the jacobian is at most 'artol'.");
    else if (converged == 2)
        sprintf(msg, "cost (one half the sum of squares of the function) is at most 'Cgoal'.");
    else if (converged == 3)
        sprintf(msg, "norm of Cost gradient is at most 'gtol'.");
    else if (converged == 4)
        sprintf(msg, "parameters changed by less than 'xtol'.");
    else if (converged == 5)
        sprintf(msg, "relative change in each of the parameters is at most 'xrtol'.");
    else if (converged == 6)
        sprintf(msg, "Cost failed to decrease by more than 'ftol' for 3 consecutive iterations.");
    else if (converged == 7)
        sprintf(msg, "relative decrease in Cost was less than 'frtol' for 3 consecutive iterations.");
    else if (info == -1)
      sprintf(msg, "Number of iterations has reached `maxiter' == %d.", nit);
	else if (info == -2)
		sprintf(msg, "Number of iterations has reached `maxfev' == %d.", nfev);
	else if (info == -3)
		sprintf(msg, "Number of iterations has reached `maxjaev' == %d.", njev);
	else if (info == -4)
		sprintf(msg, "Number of iterations has reached `maxaev' == %d.", naev);
	else if (info == -10)
		sprintf(msg, "User terminated execution via 'callback'");
	else if (info == -11)
		sprintf(msg, "NaNs produced during initial function evaluation, or subsequent Jacobian evaluations.");
    else if (info == 0)
        sprintf(msg, "Improper input parameters.");
    return msg;
}
