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
	SEXP control, SEXP rho)
{
	int		i, j;
	int		m, n, ldfjac, ldfacc;	// m - # data points, n - # pars, # ldfjac - nrow(jac)
	int		info, nfev, njev, naev;	// keep track of iterations
	int		mode;

	int		print_level, print_unit, imethod, iaccel, ibold, ibroyden, center_diff;
	int		analytic_jac, analytic_avv, damp_mode;
	double	accept, reject, lam, minlam, maxlam, avmax;

	/* Copies of inputs are stored in the C structure OS, the following are pointers to specific locations
	in OS */
	double  *par, *v, *a, *acc, *fvec, *fvec_new, **fjac, *facc, 
		*hess, *perm, *perm_t, *r, *r2, *r2_x_perm_t;
	int     *ipvt;

	// values to be returned to R.
	SEXP    eval_test;
	SEXP    sexp_diag, sexp_hess, sexp_fvec, sexp_converged, sexp_info, sexp_niter,
		sexp_message, sexp_rsstrace;
	SEXP    out, out_names;

	char    lmfun_name[8], message[256];	// stores which of the two f77 calls was used.

	int     maxfev, maxjev, maxaev;		// termination criteria

	OS = (OptStruct) R_alloc(1, sizeof(opt_struct));  // Allocate memory for structure containing inputs

	PROTECT(OS->par = duplicate(par_arg));			// Store input values in OS
	PROTECT(OS->lower = duplicate(lower_arg));
	PROTECT(OS->upper = duplicate(upper_arg));
	n = length(OS->par);							// dimension of parameter space.

	PROTECT_INDEX ipx;			// index stores the location of the "diag" element within 'control'
								// Note: Unlike 
	/* Assign the inputs into the structure, OS. Assignment procedes depending on whether variable 
	is "numeric" or a "list" */
	switch (TYPEOF(OS->par)) {
	case REALSXP:	// If par is of R-type "numeric" then do nothing.
		break;
	case VECSXP:	// If par is a list then use the list accessor/assignent functions to assign into OS, Recall that R lists are VECTORS internally.
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

	if (!isEnvironment(rho)) error("rho is not an environment!");
	OS->env = rho;	// rho is created in R by a call to new_env(), essentially blank environment.

	// checks that fcall returns a numeric vector of non-zero length, when evaluated in a blank environment
	PROTECT(eval_test = eval(OS->fcall, OS->env));
	if (!IS_NUMERIC(eval_test) || length(eval_test) == 0)
		error("evaluation of fn function returns non-sensible value!");
	m = length(eval_test);		// number of data points
	UNPROTECT(1);				// unprotect eval_test

	if (isFunction(jac)){
		PROTECT(eval_test = eval(OS->jcall, OS->env));
		if (!IS_NUMERIC(eval_test) || length(eval_test) == 0)
			error("evaluation of jac function returns non-sensible value!");
	}

	if (isFunction(acc_fun)){
		PROTECT(eval_test = eval(OS->acall, OS->env));
		if (!IS_NUMERIC(eval_test) || length(eval_test) == 0)
			error("evaluation of acc_fun function returns non-sensible value!");
	}

	ldfjac = m;		// jacobian has nrow == 'number of data points'

	// allocate memory
	a			= real_vector(n);
	v			= real_vector(n);
	par			= real_vector(n);
	acc			= real_vector(m);
	fvec		= real_vector(m);
	fvec_new	= real_vector(m);
	facc		= real_vector(m);
	fjac		= real_matrix(ldfjac, n);
	// perm		= real_vector(n * n);
	// perm_t		= real_vector(n * n);
	r			= real_vector(n * n);
	// r2			= real_vector(n * n);
	// r2_x_perm_t = real_vector(n * n);
	hess		= real_vector(n * n);

	// assign control parameters and return values to OS
	OS->ftol	= NUMERIC_VALUE(getListElement(control, "ftol"));
	OS->frtol	= NUMERIC_VALUE(getListElement(control, "frtol"));
	OS->ptol	= NUMERIC_VALUE(getListElement(control, "ptol"));
	OS->prtol	= NUMERIC_VALUE(getListElement(control, "ptol"));
	OS->gtol	= NUMERIC_VALUE(getListElement(control, "gtol"));
	OS->Cgoal	= NUMERIC_VALUE(getListElement(control, "Cgoal"));
	OS->artol   = NUMERIC_VALUE(getListElement(control, "artol"));
	OS->epsfcn1 = NUMERIC_VALUE(getListElement(control, "epsfcn1"));
	OS->epsfcn2 = NUMERIC_VALUE(getListElement(control, "epsfcn2"));
	OS->factor  = NUMERIC_VALUE(getListElement(control, "factor"));
	OS->diag	= real_matrix(n, n);


	/* Apparently, PROTECT_WITH_INDEX stores the location of the protected object so that the values
	can be overwritten when necessary. getListElement is defined in get_element.c, the function retrieves
	elements from a list by name (obviously the list must be named for this to work). Internally Matrices
	are of type REALSXP (probably) and so all that needs to be done is to write a conversion between
	the array passed too and from FORTRAN and the SEXP class object stored in OS->diag */
	PROTECT_WITH_INDEX(sexp_diag = getListElement(control, "diag"), &ipx);
	switch (TYPEOF(sexp_diag)) {
	case REALSXP:
		// having created a copy of "diag" assign pointers to each element of diag into OS.
		if (length(sexp_diag) == n * n) {
			REPROTECT(sexp_diag = duplicate(sexp_diag), ipx);
			for (i = 0; i < n; i++) {
				for (j = 0; j < n; j++) {
					OS->diag[i][j] = NUMERIC_POINTER(sexp_diag)[i + j * n];
				}
			}
			mode = 2;	// What's this for? (numeric code for outcome of length(sexp_diag) == n?)
		}
		else {
			REPROTECT(sexp_diag = NEW_NUMERIC(n * n), ipx);
			mode = 1;
		}

		sexp_diag = SET_ARRAY_DIM(sexp_diag, n, n);		// set dim attribute of output array (2016-03-30-01:45)
		break;
	case VECSXP:
		#ifndef NO_VECSXP
		// Not exactly sure how the above approach translates to the VECSXP case
		if (length(sexp_diag) == n) {
			REPROTECT(sexp_diag = duplicate(sexp_diag), ipx);
			for (i = 0; i < n; i++) {
				SET_VECTOR_ELT(sexp_diag, i, AS_NUMERIC(VECTOR_ELT(sexp_diag, i)));
				OS->diag[i] = NUMERIC_VALUE(VECTOR_ELT(sexp_diag, i));
			}
			mode = 2;
		}
		else {
			REPROTECT(sexp_diag = NEW_LIST(n), ipx);
			for (i = 0; i < n; i++)
				SET_VECTOR_ELT(sexp_diag, i, NEW_NUMERIC(1));
			mode = 1;
		}
		#else 
		error("control$diag must be of class 'numeric'!");
		#endif
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
	OS->nprint = INTEGER_VALUE(getListElement(control, "nprint"));
	if (OS->nprint > 0)
		OS->nprint = 1;
	/* The memory location reffered to by a pointer is the first "block" of memory allocated to that
	variable. The references declared at the top of the function definition are used to refer to the
	last "block" of memory allocated to each variable. */
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

	// initialise number of iterations and iflag.
	OS->niter	  = 0;
	OS->converged = 0;

	// Set flags that dictate the method
	info		 = 0;
	print_level  = 4;
	print_unit	 = 2;
	imethod		 = 1;
	iaccel		 = 0;
	ibold		 = 0;
	ibroyden	 = 0;
	center_diff  = 0;
	analytic_jac = 0;
	analytic_avv = 0;
	damp_mode	 = 1;
	accept		 = 2.0;
	reject		 = 5.0;
	maxlam		 = 100;
	minlam		 = 0.1;
	avmax		 = 100;


	Rprintf("Calling geodesiclm!\n");
	/*  */
	F77_CALL(geodesiclm)(&fcn_lm, &fcn_ja, &fcn_ac, par, fvec, fjac, &n, &m,
		&fcn_call, &info, &analytic_jac, &analytic_avv, &center_diff, &OS->epsfcn1, &OS->epsfcn2,
		OS->diag, &damp_mode, &OS->niter, &nfev, &njev, &naev, &OS->maxiter, &maxfev, &maxjev, &maxaev,
		&maxlam, &minlam, &OS->artol, &OS->Cgoal, &OS->gtol, &OS->ptol, &OS->prtol, &OS->ftol, &OS->frtol,
		&OS->converged, &print_level, &print_unit, &imethod, &iaccel, &ibold, &ibroyden,
		&OS->factor, &accept, &reject, &avmax);

	UNPROTECT(1);
	strcpy(lmfun_name, "geodesiclm");

	// Store diagnostic message for regression output
	fcn_message(message, OS->converged, info, n, OS->niter, nfev, njev, naev);
	if ( (OS->converged < 1 || 8 < OS->converged) && (info < -12 || info > 1) )
		warning("%s: info = %d. %s\n\n", lmfun_name, info, message);

	Rprintf("Final info = %i \n", info);
	Rprintf("Final convergence = %i \n", OS->converged);
	Rprintf("Hessian matrix calculations \n");
	// allocate memory to store the hessian matrix, hessian is stored in a compressed format.
	PROTECT(sexp_hess = NEW_NUMERIC(n*n));
	for (j = 0; j < n; j++)
		for (i = 0; i < n; i++) {
			r[j*n + i] = (i <= j) ? fjac[i][j] : 0;
			Rprintf("r[%i*n + %i] = %g \n", j, i, r[j*n + i]);
		}

	/*  perm %*% t(r) %*% r %*% t(perm)  *
	*    |       |___r2__|         |    *
	*    |           |_r2_x_perm_t_|    *
	*    |_______hess_______|           */
	/*
	Rprintf("transpose \n");
	transpose(perm, n, n, perm_t);
	Rprintf("crossprod \n");
	crossprod(r, n, n, r, n, n, r2);
	Rprintf("matprod1 \n");
	matprod(r2, n, n, perm_t, n, n, r2_x_perm_t);
	Rprintf("matprod2 \n");
	matprod(perm, n, n, r2_x_perm_t, n, n, hess);
	*/

	/*	t(r) %*% r 
		|__hess__|	*/
	Rprintf("hessian");
	crossprod(r, n, n, r, n, n, hess);

	Rprintf("Assign to sexp_hess \n");
	for (i = 0; i < n*n; i++)
		NUMERIC_POINTER(sexp_hess)[i] = hess[i];

	Rprintf("Assign to sexp_fvec \n");
	// initialise residual vector and create reference to the end of the vector.
	PROTECT(sexp_fvec = NEW_NUMERIC(m));
	for (i = 0; i < m; i++)
		NUMERIC_POINTER(sexp_fvec)[i] = fvec[i];

	// initialise vector of the rss at each iteration and create reference to the end of the vector.
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
				NUMERIC_POINTER(sexp_diag)[i + i * n]= OS->diag[i][i];
		}
	}
	else {
		#ifndef NO_VECSXP
		for (i = 0; i < n; i++)
			NUMERIC_POINTER(VECTOR_ELT(sexp_diag, i))[0] = OS->diag[i];
		#else 
		error("control$diag must be of class 'numeric'!");
		#endif
	}

	Rprintf("Creating output list\n");
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
	UNPROTECT(13);
	free(fjac);
	free(OS->diag);

	// aaaand done
	return out;
}