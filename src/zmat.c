/*
 *  zmat.c: Sample just Zs
 *  Copyright (C) 2016 Charles C. Berry <cberry@ucsd.edu>
 *
 *  This program is free software; you can redistribute it
 *  and/or modify it under the terms of the GNU General
 *  Public License as published by the Free Software
 *  Foundation; either version 2 of the License, or (at your
 *  option) any later version.
 *
 *  This program is distributed in the hope that it will be
 *  useful, but WITHOUT ANY WARRANTY; without even the
 *  implied warranty of MERCHANTABILITY or FITNESS FOR A
 *  PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 *
 *  You should have received a copy of the GNU General
 *  Public License along with this program; if not, a copy
 *  is available at http://www.r-project.org/Licenses/
 *  
 *  This provides the function innerLoop and support
 *  functions for it.  It is intended that this functon be
 *  called by R.  See the R code and help pages for details
 *  on setting up the objects used.
 */
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

void zmat(
	  double *prw,
	  double *pz,
	  double *eps,
	  double *V,
	  double *eta,
	  double *omdp,
	  int *w,
	  int *wp,
	  int *n,
	  int *T,
	  int *ka,
	  int *ko,
	  int *ndat,
	  int *zy,
	  double *etaomdp,
	  double *workT
	  ){

  /* calc etaomdp and prw */
  matmult(eta, *T, *ka, omdp, *ka, *ko, etaomdp);

  /* rowsums */
  for (int t=0;t<*T;t++){
    double x=0.0;
    for (int k=0;k<*ka;k++) x+=etaomdp[t+k**T];
    prw[t]=x;
  }

  /* probzv */

  probzv(V, pz, T);

  /* zy */

  zysum(prw, pz, eps, eta, omdp, w, wp, n, T, ka, ko, ndat, zy, etaomdp, workT);

}
