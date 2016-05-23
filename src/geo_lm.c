#ifdef LOCAL_HEADERS
#include "R.h"
#include "Rdefines.h"
#include "GeodesicLM.h"

#else
#include <R.h>
#include <Rdefines.h>
#include "GeodesicLM.h"

#endif

OptStruct OS;

SEXP geo_lm(SEXP par_arg, SEXP lower_arg, SEXP upper_arg, SEXP fn, SEXP jac, SEXP acc_fun,
	SEXP callback, SEXP control, SEXP rho)
{
	/*************************************************************************************/
	// Variable Declarations and some memory allocation
	
	int		i, j;
	int		m, n, ldfjac, ldfacc;	// m - # data points, n - # pars, # ldfjac - nrow(jac)
	int		info, nfev, njev, naev;	// keep track of iterations and evaluations
	int		mode;

	int		print_level, print_unit, imethod, iaccel, ibold, ibroyden, center_diff;
	int		analytic_jac, analytic_avv, damp_mode;
	double	accept, reject, lam, minlam, maxlam, avmax;

	/* inputs and outputs are stored in the structure pointed to by 'OS', seperate instances of
	the parameter vector, jacobian, hessian etc. are passed to the geodesiclm subroutine. */
	double  *par, *v, *a, *acc, *fvec, *fvec_new, *fjac, *facc, *hess, *r;
	int     *ipvt;

	// values to be returned to R.
	SEXP    eval_test;
	SEXP    sexp_diag, sexp_hess, sexp_fvec, sexp_converged, sexp_info, sexp_niter,
		sexp_message, sexp_rsstrace;
	SEXP    out, out_names;

	char    lmfun_name[8], message[256];	// left-over from minpack.lm code

	int     maxfev, maxjev, maxaev;		// termination criteria

	// Allocate memory for and assign inputs to opt_struct pointed to by OS.
	OS = (OptStruct) R_alloc(1, sizeof(opt_struct)); 
	PROTECT(OS->par = duplicate(par_arg));
	PROTECT(OS->lower = duplicate(lower_arg));
	PROTECT(OS->upper = duplicate(upper_arg));
	n = length(OS->par);						// dimension of parameter space.

	PROTECT_INDEX ipx;			// index stores the location of the "diag" element within 'control'


	/******************************************************************************/
	// Input validation 

	/* The following code appears to be a safe-guard against the user supplying a list containing 
	non-numeric elements, or elements that cannot be coerced to numeric. */
	switch (TYPEOF(OS->par)) {
	case REALSXP:
		break;
	case VECSXP:
		for (i = 0; i < n; i++)
			SET_VECTOR_ELT(OS->par, i, AS_NUMERIC(VECTOR_ELT(OS->par, i)));
		break;
	default:
		error("`par' that you provided is non-list and non-numeric!");
	}

	switch (TYPEOF(OS->lower)) {
	case REALSXP:
		break;
	case VECSXP:
		for (i = 0; i < n; i++)
			SET_VECTOR_ELT(OS->lower, i, AS_NUMERIC(VECTOR_ELT(OS->lower, i)));
		break;
	default:
		error("`lower' that you provided is non-list and non-numeric!");
	}

	switch (TYPEOF(OS->upper)) {
	case REALSXP:
		break;
	case VECSXP:
		for (i = 0; i < n; i++)
			SET_VECTOR_ELT(OS->upper, i, AS_NUMERIC(VECTOR_ELT(OS->upper, i)));
		break;
	default:
		error("`upper' that you provided is non-list and non-numeric!");
	}

	if (!isFunction(fn)) error("fn is not a function!");
	PROTECT(OS->fcall = lang2(fn, OS->par));	// assign lhs function into OS (guessing that's what lang2 does)
	
	if (!isNull(jac))
		PROTECT(OS->jcall = lang2(jac, OS->par));
	if (!isNull(acc_fun))
		PROTECT(OS->acall = lang2(acc_fun, OS->par));
	if (!isNull(callback))
		PROTECT(OS->callback = lang2(callback, OS->par));

	if (!isEnvironment(rho)) error("rho is not an environment!");
	OS->env = rho;	// rho is created in R by a call to new_env(), essentially blank environment.

	// checks that fcall returns a numeric vector of non-zero length, when evaluated in a blank environment
	PROTECT(eval_test = eval(OS->fcall, OS->env));
	if (!IS_NUMERIC(eval_test) || length(eval_test) == 0)
		error("evaluation of fn function returns non-sensible value!");
	m = length(eval_test);		// number of data points
	UNPROTECT(1);				// unprotect eval_test


	/************************************************************************************/
	// Initialisation and more memory allocation

	/* set flags to calculate jacobian and acceleration using finite differences, unless functions are
	supplied by the user. */
	analytic_jac = 0;
	analytic_avv = 0;

	if (isFunction(jac)){
		analytic_jac = 1;
		PROTECT(eval_test = eval(OS->jcall, OS->env));
		if (!IS_NUMERIC(eval_test) || length(eval_test) == 0)
			error("evaluation of jac function returns non-sensible value!");
	}
	UNPROTECT(1);				// unprotect eval_test

	if (isFunction(acc_fun)){
		analytic_avv = 1;
		PROTECT(eval_test = eval(OS->acall, OS->env));
		if (!IS_NUMERIC(eval_test) || length(eval_test) == 0)
			error("evaluation of acc_fun function returns non-sensible value!");
	}
	UNPROTECT(1);				// unprotect eval_test

	if (isFunction(callback)) {
		PROTECT(eval_test = eval(OS->callback, OS->env));
		if (!IS_NUMERIC(eval_test) || length(eval_test) == 0)
			error("evaluation of callback function returns non-sensible value!");
	}
	UNPROTECT(1);				// unprotect eval_test

	ldfjac = m;		// When is ldfjac not equal to 'm'?

	// allocate memory for variables to be passed to geodesiclm subroutine
	a			= real_vector(n);
	v			= real_vector(n);
	par			= real_vector(n);
	acc			= real_vector(m);
	fvec		= real_vector(m);
	fvec_new	= real_vector(m);
	facc		= real_vector(m);
	fjac		= real_vector(ldfjac * n);
	r			= real_vector(n * n);
	hess		= real_vector(n * n);

	// read parameters from `control` and return values to seperate variables in OS
	// getListElement:		retrieves elements from a list by name (obviously the list must
	//						be named for this to work)
	OS->ftol			= NUMERIC_VALUE(getListElement(control, "ftol"));
	OS->frtol			= NUMERIC_VALUE(getListElement(control, "frtol"));
	OS->xtol			= NUMERIC_VALUE(getListElement(control, "xtol"));
	OS->xrtol			= NUMERIC_VALUE(getListElement(control, "xrtol"));
	OS->gtol			= NUMERIC_VALUE(getListElement(control, "gtol"));
	OS->Cgoal			= NUMERIC_VALUE(getListElement(control, "Cgoal"));
	OS->artol			= NUMERIC_VALUE(getListElement(control, "artol"));
	OS->epsfcn1			= NUMERIC_VALUE(getListElement(control, "epsfcn1"));
	OS->epsfcn2			= NUMERIC_VALUE(getListElement(control, "epsfcn2"));
	OS->damp_mode		= INTEGER_VALUE(getListElement(control, "damp_mode"));
	OS->center_diff		= INTEGER_VALUE(getListElement(control, "center_diff")); //coercion to int from LGLSXP?
	OS->diag			= real_vector(n * n);

	// PROTECT_WITH_INDEX:	Stores the location of the protected object so that the values
	//						can be overwritten when necessary. (not really sure how/why this works)
	// REPROTECT:			Not sure how this works or why it's used.
	// duplicate control$diag to output and make a local copy to pass to geodesiclm subroutine
	PROTECT_WITH_INDEX(sexp_diag = getListElement(control, "diag"), &ipx);
	switch (TYPEOF(sexp_diag)) {
	case REALSXP:
		if (length(sexp_diag) == n * n) {
			REPROTECT(sexp_diag = duplicate(sexp_diag), ipx);
			for (i = 0; i < n; i++) {
				for (j = 0; j < n; j++) {
					OS->diag[i + n * j] = NUMERIC_POINTER(sexp_diag)[i + j * n];
				}
			}
			mode = 2;	// numeric code for outcome of length(sexp_diag) == n*n?
		}
		else {
			REPROTECT(sexp_diag = NEW_NUMERIC(n * n), ipx);
			mode = 1;
		}
		break;
	case VECSXP:
		if (length(sexp_diag) == n * n) {
			REPROTECT(sexp_diag = duplicate(sexp_diag), ipx);
			for (i = 0; i < n; i++) {
				for (j = 0; j < n; j++) {
					SET_VECTOR_ELT(sexp_diag, i + j*n, AS_NUMERIC(VECTOR_ELT(sexp_diag, i + j*n)));
					OS->diag[i + j*n] = NUMERIC_VALUE(VECTOR_ELT(sexp_diag, i + j*n));
				}
			}
			mode = 2;
		}
		else {
			REPROTECT(sexp_diag = NEW_LIST(n*n), ipx);
			for (i = 0; i < n; i++) {
				for (j = 0; j < n; j++) {
					SET_VECTOR_ELT(sexp_diag, i + j*n, NEW_NUMERIC(1));
				}
			}
			mode = 1;
		}
		break;
	default:
		error("`diag' that you provided is non-list and non-numeric!");
	}

	// assign more control parameters.
	maxfev = INTEGER_VALUE(getListElement(control, "maxfev"));
	maxjev = INTEGER_VALUE(getListElement(control, "maxjev"));
	maxaev = INTEGER_VALUE(getListElement(control, "maxaev"));
	OS->maxiter = INTEGER_VALUE(getListElement(control, "maxiter"));
	if (OS->maxiter > 1024) {
		OS->maxiter = 1024;
		warning("resetting `maxiter' to 1024!");
	}
	OS->print_level		= INTEGER_VALUE(getListElement(control, "print_level"));
	OS->imethod			= INTEGER_VALUE(getListElement(control, "imethod"));
	OS->iaccel			= INTEGER_VALUE(getListElement(control, "iaccel"));
	OS->ibold			= INTEGER_VALUE(getListElement(control, "ibold"));
	OS->ibroyden		= INTEGER_VALUE(getListElement(control, "ibroyden"));
	OS->initial_factor	= NUMERIC_VALUE(getListElement(control, "initial_factor"));
	OS->accept			= NUMERIC_VALUE(getListElement(control, "accept"));
	OS->reject			= NUMERIC_VALUE(getListElement(control, "reject"));
	OS->avmax			= NUMERIC_VALUE(getListElement(control, "avmax"));
	
	if (IS_NUMERIC(OS->par)) {
		for (i = 0; i < n; i++) // Set 'i' to 'n-1' essentially?
			if (R_FINITE(NUMERIC_POINTER(OS->par)[i]))
				par[i] = NUMERIC_POINTER(OS->par)[i];
			else
				error("Non-finite (or null) value for a parameter specified!");
	}
	else {
		for (i = 0; i < n; i++)
			if (R_FINITE(NUMERIC_VALUE(VECTOR_ELT(OS->par, i))))
				par[i] = NUMERIC_VALUE(VECTOR_ELT(OS->par, i));
			else
				error("Non-finite (or null) value for a parameter specified!");
	}

	// initialise number of iterations and convergence flag.
	OS->niter	  = 0;
	OS->converged = 0;

	// Set flags that dictate the method
	info		 = 0;		// User termination flag
	print_unit	 = 6;		// prints messages from geodesiclm to the console
	maxlam		 = -100;	// Negative values provide no limit to the value of lambda
	minlam		 = -0.1;

	// geodesiclm performs the optimisation.
	F77_CALL(geodesiclm)(&fcn_lm, &fcn_ja, &fcn_ac, par, fvec, fjac, &n, &m,
		&fcn_call, &info, &analytic_jac, &analytic_avv, &OS->center_diff, &OS->epsfcn1, &OS->epsfcn2,
		OS->diag, &OS->damp_mode, &OS->niter, &nfev, &njev, &naev, &OS->maxiter, &maxfev, &maxjev, &maxaev,
		&maxlam, &minlam, &OS->artol, &OS->Cgoal, &OS->gtol, &OS->xtol, &OS->xrtol, &OS->ftol, &OS->frtol,
		&OS->converged, &OS->print_level, &print_unit, &OS->imethod, &OS->iaccel, &OS->ibold, &OS->ibroyden,
		&OS->initial_factor, &OS->accept, &OS->reject, &OS->avmax);

	strcpy(lmfun_name, "geodesiclm");

	// Store diagnostic message for regression output, based on values of convergence/stopping flags
	fcn_message(message, OS->converged, info, n, OS->niter, nfev, njev, naev);
	if ( (OS->converged < 1 || 8 < OS->converged) && (info < -12 || info > 1) )
		warning("%s: info = %d. %s\n\n", lmfun_name, info, message);

	// Calculations for the hessian matrix (isn't this calculated in geodesiclm already?).
	PROTECT(sexp_hess = NEW_NUMERIC(n*n));
	for (j = 0; j < n; j++){
		for (i = 0; i < n; i++) {
			r[j*n + i] = (i <= j) ? fjac[i + ldfjac*j] : 0;
		}
	}

	/*			 t(r) %*% r				   *
	*    |      |___hess___|         |    */
	crossprod(r, n, n, r, n, n, hess);

	for (i = 0; i < n*n; i++)
		NUMERIC_POINTER(sexp_hess)[i] = hess[i];

	// Copy values of fvec to output REALSXP.
	PROTECT(sexp_fvec = NEW_NUMERIC(m));
	for (i = 0; i < m; i++)
		NUMERIC_POINTER(sexp_fvec)[i] = fvec[i];

	// Copy values of rsstrace to output REALSXP.
	PROTECT(sexp_rsstrace = NEW_NUMERIC(OS->niter));
	for (i = 0; i < OS->niter; i++)
		NUMERIC_POINTER(sexp_rsstrace)[i] = OS->rsstrace[i];

	PROTECT(sexp_converged = NEW_INTEGER(1));
	INTEGER_POINTER(sexp_converged)[0] = OS->converged;

	PROTECT(sexp_info = NEW_INTEGER(1));
	INTEGER_POINTER(sexp_info)[0] = info;

	PROTECT(sexp_niter = NEW_INTEGER(1));
	INTEGER_POINTER(sexp_niter)[0] = OS->niter - 1;

	PROTECT(sexp_message = NEW_STRING(1));
	SET_STRING_ELT(sexp_message, 0, mkChar(message));
	if (IS_NUMERIC(sexp_diag)) {
		for (i = 0; i < n; i++) {
				NUMERIC_POINTER(sexp_diag)[i + i * n]= OS->diag[i + n*i];
		}
	}
	else {
		for (i = 0; i < n; i++)
			NUMERIC_POINTER(VECTOR_ELT(sexp_diag, i))[0] = OS->diag[i];
	}


	// create a list storing the output to be returned to R.
	PROTECT(out = NEW_LIST(9));
	SET_VECTOR_ELT(out, 0, OS->par);
	SET_VECTOR_ELT(out, 1, sexp_hess);
	SET_VECTOR_ELT(out, 2, sexp_fvec);
	SET_VECTOR_ELT(out, 3, sexp_converged);
	SET_VECTOR_ELT(out, 4, sexp_info);
	SET_VECTOR_ELT(out, 5, sexp_message);
	SET_VECTOR_ELT(out, 6, sexp_diag);
	SET_VECTOR_ELT(out, 7, sexp_niter);
	SET_VECTOR_ELT(out, 8, sexp_rsstrace);

	// assign names to each element of the output
	PROTECT(out_names = NEW_STRING(9));
	SET_STRING_ELT(out_names, 0, mkChar("par"));
	SET_STRING_ELT(out_names, 1, mkChar("hessian"));
	SET_STRING_ELT(out_names, 2, mkChar("fvec"));
	SET_STRING_ELT(out_names, 3, mkChar("converged"));
	SET_STRING_ELT(out_names, 4, mkChar("info"));
	SET_STRING_ELT(out_names, 5, mkChar("message"));
	SET_STRING_ELT(out_names, 6, mkChar("diag"));
	SET_STRING_ELT(out_names, 7, mkChar("niter"));
	SET_STRING_ELT(out_names, 8, mkChar("rsstrace"));

	SET_NAMES(out, out_names);

	// remove protections
	UNPROTECT(18);

	// aaaand done
	return out;
}