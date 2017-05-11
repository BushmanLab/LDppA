#include <R.h>
#include <Rmath.h>
#include <R_ext/Random.h>
#include <R_ext/Print.h>
#include <R_ext/BLAS.h>

void matmult(double *x, int nrx, int ncx,
	     double *y, int nry, int ncy, double *z);
void probzv( double *V, double *Z, int *T);
void zysum(
	   double *prw, double *pz, double *eps, double *eta,
	   double *omdp, int *w, int *wp, int *n,
	   int *T, int *ka, int *ko, int *ndat, int *zy,
	   double *etaomdp, double *workT);
