/* 
 * utils.c: Support funs for ectc
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
#include "ectc.h"
/* utils */

void matmult(double *x, int nrx, int ncx,
	     double *y, int nry, int ncy, double *z)
{
  char *transa = "N", *transb = "N";
  double one = 1.0, zero = 0.0;

  F77_CALL(dgemm)(transa, transb, &nrx, &ncy, &ncx, &one,
		  x, &nrx, y, &nry, &zero, z, &nrx);
}

// //* As per Equation \ref{eq:prz}
void probzv( double *V, double *Z, int *T){
  double usable=1.0;
  for (int t = 0; t < *T; t++) {
    Z[t]=V[t]*usable;
    usable*=fmax2(1-V[t],0.0);
  }
}
void zysum(
	   double *prw,
	   double *pz,
	   double *eps,
	   double *eta,
	   double *omdp,
	   int *w,
	   int *wp,
	   int *n,
	   int *T,
	   int *ka,
	   int *ko,
	   int *ndat,
	   /* result */
	   int *zy,
	   /* workspace */
	   double *etaomdp,
	   double *workT
	   ){
  // //* Refer to Section \ref{sec:zdist}
  for (int idat=0;idat<*ndat;idat++){
    double biglog=R_NegInf;
    for (int t=0;t<*T;t++){
      int one=1L;
      // //* Factor Equation \ref{eq:prw.z5} as multinomial times negative binomial
      // Negative Binomial Part
      // Note: eps==0 is limiting case
      // Note: prob parm needs to be 1-nbparm given R convention
      double x = (eps[0]==0) ? 0.0 :
	dnbinom(wp[idat],eps[0],1-prw[t]/(eps[1]+prw[t]),one);
      for (int k=0;k<*ko;k++) x+= log(etaomdp[t+k**T]/prw[t]) *
				w[idat+k**ndat];
      if (biglog<x) biglog=x;
      workT[t]=x;
    }
    double prTot=0.0;
    for (int t=0;t<*T;t++){
      workT[t]=exp(workT[t]-biglog)*pz[t];
      prTot+=workT[t];
    }
    // //* Multinomial part of Equation \ref{eq:prw.z5}:
    for (int t=0;t<*T;t++) workT[t]/=prTot;
    rmultinom((int) n[idat],workT,(int) T[0], zy+idat**T);
  }  
}
// sample pz from dirichlet


void samplePz(
	      int *T,            // number of compositions
	      int *ndat,         // number of data patterns
	      double *alpha,     // prior for V
	      int *zy,           // T by ndat table of assignments
	      double *pz,        // prob Z = t
	      double *workT     // workspace
	      ){
  // accum zsums
	      
  for (int i=0; i<*T; i++) workT[i]=*alpha;
	      
  for (int j=0; j<*ndat; j++)
    for (int i=0; i<*T; i++) workT[i]+=(double) zy[i+*T*j];
	      
  double gsum=0.0;
  for (int i=0; i<*T; i++){
    pz[i] = rgamma(workT[i], 1.0);
    gsum+=pz[i];
  }
	      
  for (int i=0; i<*T; i++) pz[i]/=gsum;
	      
	      
}
