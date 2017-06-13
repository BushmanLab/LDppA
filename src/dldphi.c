#include <R.h>
#include <Rinternals.h>

SEXP dldphi( SEXP phi, SEXP w, SEXP op ){
  int lenphi = length(phi);
  int neta = lenphi+1;

  double etadenom = 1.0;
  for (int i = 0 ; i< lenphi; i++) etadenom += exp(REAL(phi)[i]);

  double *eta = Calloc(neta, double);

  eta[0]=1.0/etadenom;
  for (int i = 0; i<length(phi); i++) eta[i+1]=exp(REAL(phi)[i])/etadenom;

  double *eodp = Calloc( neta, double); // later becomes p
  double *rowsm = Calloc( neta, double);
  double pdenom = 0.0;


  for (int i = 0; i<neta; i++){
    double rowsmi = 0.0;
    for (int j = 0;j<neta; j++) rowsmi += REAL(op)[neta*j + i ];
    rowsm[i]=rowsmi;
  }



  for (int j = 0;j<neta; j++){  
    double eodpj = 0.0;
    for (int i = 0; i<neta; i++) eodpj += eta[i] * REAL(op)[neta*j + i ];
    pdenom+=eodpj;
    eodp[j]=eodpj;
  }

  for (int i = 0; i<neta; i++) eodp[i]/=pdenom; // now eodp is p

  SEXP dlkdphi = PROTECT(allocVector(REALSXP, lenphi));

  for (int k = 0; k < lenphi; k++){
    double res = 0.0;
    for (int i = 0; i<neta; i++){
      for (int j = 0; j < neta; j++){
	res += (REAL(w)[i]/eodp[i]) * // dllk.dp
	  (REAL(op)[ i * neta + j]/pdenom - eodp[i] * rowsm[j] / pdenom ) * // dp.deta 
	  ( (j == k +1 ? eta[j] : 0.0 ) + eta[j]*exp(REAL(phi)[k]/etadenom));
      }
    }
    REAL(dlkdphi)[k] = res;
  }

  Free(eta); Free(eodp); Free(rowsm);
  UNPROTECT(1);

  return dlkdphi;

}
