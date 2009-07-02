/* Define some global variables and functions that operate on some of them
Similar to deSolve_utils.c from package deSolve*/

#include <R.h>
#include <Rdefines.h>
#include "steady.h"

/*==================================================
some functions for keeping track of how many SEXPs
are PROTECTed, and UNPROTECTing them in the
case of a fortran stop.
==================================================*/
long int N_Protected;

void init_N_Protect(void) { N_Protected = 0; }

void incr_N_Protect(void) { N_Protected++; }

void unprotect_all(void) { UNPROTECT((int) N_Protected); }

void my_unprotect(int n)
{
    UNPROTECT(n);
    N_Protected -= n;
}

/* Globals :*/ 

SEXP steady_deriv_func;
SEXP steady_jac_func;
SEXP steady_envir;
SEXP st_gparms;

SEXP stsparse_deriv_func;
SEXP stsparse_jac_func;
SEXP stsparse_envir;

/* runsteady */
SEXP lsode_deriv_func;
SEXP lsode_jac_func;
SEXP lsode_envir;

/*==================================================
Parameter initialisation functions
===================================================*/

void initParms(SEXP Initfunc, SEXP Parms) {

  if (inherits(Initfunc, "NativeSymbol"))  {
    init_func *initializer;

    PROTECT(st_gparms = Parms);     incr_N_Protect();
    initializer = (init_func *) R_ExternalPtrAddr(Initfunc);
    initializer(Initstparms);
  }

}

void Initstparms(int *N, double *parms)
{
  int i, Nparms;

  Nparms = LENGTH(st_gparms);
  if ((*N) != Nparms)
    {
      PROBLEM "Confusion over the length of parms"
      ERROR;
    } 
  else
    {
      for (i = 0; i < *N; i++) parms[i] = REAL(st_gparms)[i];
    }
}
  
/*==================================================
 output initialisation function

 out and ipar are used to pass output variables
 (number set by nout) followed by other input
 by R-arguments rpar, ipar
 ipar[0]: number of output variables,
 ipar[1]: length of rpar,
 ipar[2]: length of ipar
===================================================*/

void initOut(int isDll, int neq, SEXP nOut, SEXP Rpar, SEXP Ipar) {

  int j;
  nout   = INTEGER(nOut)[0];    /* number of output variables */
  if (isDll)                    /* function is a dll */
  {
   if (nout > 0) isOut = 1;
   ntot  = neq + nout;          /* length of yout */
   lrpar = nout + LENGTH(Rpar); /* length of rpar; LENGTH(Rpar) is always >0 */
   lipar = 3 + LENGTH(Ipar);    /* length of ipar */

  } else                              /* function is not a dll */
  {
   isOut = 0;
   ntot = neq;
   lipar = 1;
   lrpar = 1;
  }

   out   = (double *) R_alloc(lrpar, sizeof(double));
   ipar  = (int *)    R_alloc(lipar, sizeof(int));

   if (isDll ==1)
   {
    ipar[0] = nout;              /* first 3 elements of ipar are special */
    ipar[1] = lrpar;
    ipar[2] = lipar;
    /* other elements of ipar are set in R-function via argument *ipar* */
    for (j = 0; j < LENGTH(Ipar);j++) ipar[j+3] = INTEGER(Ipar)[j];

    /* first nout elements of rpar reserved for output variables
      other elements are set in R-function lsodx via argument *rpar* */
    for (j = 0; j < nout; j++)        out[j] = 0.;
    for (j = 0; j < LENGTH(Rpar);j++) out[nout+j] = REAL(Rpar)[j];
   }

}


