#ifndef GEODESIC_LM_HEADER_INCLUDED
#define GEODESIC_LM_HEADER_INCLUDED 1
#define NO_VECSXP		// if defined then the VECSXP case is ignored in some places

// #define LOCAL_HEADERS	// if defined then use local copies of the R header files (to appease intellisense)

#ifdef LOCAL_HEADERS
#include "R.h"
#include "Rdefines.h"
#include "GeodesicLM.h"

#else
#include <R.h>
#include <Rdefines.h>
#include "GeodesicLM.h"

#endif

typedef struct opt_struct {
	SEXP par;
	SEXP lower;
	SEXP upper;
	SEXP fcall;
	SEXP jcall;
	SEXP acall;
	SEXP env;
	double ftol;
	double frtol;
	double Cgoal;
	double artol;
	double gtol;
	double ptol;
	double prtol;
	double epsfcn1;
	double epsfcn2;
	double *diag;
	double factor;
	double accept;
	double reject;
	int nprint;
	int maxiter;
	int niter;
	int converged;
	double rsstrace[1024];
}opt_struct, *OptStruct;

void fcn_lm(int *m, int *n, double *par, double *fvec);
void fcn_ja(int *ldfjac, int *n, double *par, double *fjac);
void fcn_ac(int *ldfacc, int *n, double *par, double *v, double *facc);
void fcn_call(int *m, int *n, double *par, double *v, double *a, double *fvec, double *fjac,
	double *acc, double *lam, double *dtd, double *fvec_new, int *accepted, int *info);

void F77_NAME(geodesiclm)(void(*fcn_lm)(int *m, int *n, double *par, double *fvec),
	void(*fcn_ja)(int *ldfjac, int *n, double *par, double *fjac), 
	void(*fcn_ac)(int *ldfacc, int *n, double *par, double *v, double *facc),
	double *par, double *fvec, double **fjac, int *n, int *m, 
	void(*fcn_call)(int *m, int *n, double *par, double *v, double *a,
		double *fvec, double *fjac, double *acc,
		double *lam, double *dtd, double *fvec_new, int *accepted, int *info),
	int *info, int *analytic_jac, int *analytic_Avv, int *center_diff, double *epsfcn1, double *epsfcn2, 
	double **diag, int *damp_mode, int *niters, int *nfev, int *njev, int *naev,
	int *maxiters, int *maxfev, int *maxjev, int *maxaev, double *maxlam, double *minlam,
	double *artol, double *Cgoal, double *gtol, double *ptol, double *prtol, double *ftol, double *frtol,
	int *converged, int *print_level, int *print_unit, 
	int *imethod, int *iaccel, int *ibold, int *ibroyden, double *factor, double *accept, double *reject, 
	double *avmax);

SEXP	getListElement(SEXP list, char *str);
SEXP	SET_ARRAY_DIM(SEXP x, int m, int n);
double	*real_vector(int n);
int		*int_vector(int n);

void transpose(double *x, int nrx, int ncx, double *y);
void matprod(double *x, int nrx, int ncx, double *y, int nry, int ncy, double *z);
void crossprod(double *x, int nrx, int ncx, double *y, int nry, int ncy, double *z);
char *fcn_message(char *msg, int converged, int info, int n, int nit, int nfev, int njev, int naev);


SEXP geo_lm(SEXP par_arg, SEXP lower_arg, SEXP upper_arg,
	SEXP fn, SEXP jac, SEXP acc, SEXP control, SEXP rho);

extern OptStruct OS;

#endif /* GEODESIC_LM_HEADER_INCLUDED */
