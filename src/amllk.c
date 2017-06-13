#include <R.h>
#include <Rinternals.h>

SEXP amllk( SEXP phi, SEXP w, SEXP omegaPsi){

  SEXP res = PROTECT(allocVector(REALSXP, 1));

  int lenphi = length(phi);
  int leneta = lenphi + 1;

  double *eta, *etaOP, totexpEta = 1.0,
    maxphi = 0.0, sumet = 0.0;


  for (int i=0; i<lenphi; i++){
    double p = REAL(phi)[i];
    if (p>maxphi) maxphi = p;
  }

  eta =  Calloc( leneta , double );
  eta[0] = exp(-maxphi);
  totexpEta = eta[0];
  for (int i=0; i<lenphi; i++){
    eta[i+1] = exp(REAL(phi)[i] - maxphi); 
    totexpEta += eta[i+1];
  }

  for (int i=0; i<leneta; i++) eta[i]/=totexpEta;

  etaOP = Calloc( leneta, double);

  for (int i =0; i<leneta; i++){
    double et = 0.0;
    for (int j = 0; j < leneta; j++)
      et += eta[ j ] * REAL(omegaPsi)[ j + leneta * i];
    etaOP[i]=et;
    sumet+=et;
  }

  double rs = 0.0;
  for (int i=0; i<leneta; i++) rs+= REAL(w)[i]*log(etaOP[i]/sumet);

  Free(eta);
  Free(etaOP);

  REAL(res)[0] = rs;

  UNPROTECT(1);

  return res;

}
