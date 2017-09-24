#include <RcppArmadillo.h>
using namespace Rcpp;


RcppExport SEXP fixedEtaEMCall(
			       SEXP lkMatrix,
			       SEXP probg,
			       SEXP ritab,
			       SEXP alpha,
			       SEXP weights,
			       SEXP nreps
			       ){
  double alfa = Rcpp::as<double>(alpha);
  int nrp =  Rcpp::as<int>(nreps);
  arma::mat lkm = Rcpp::as<arma::mat>(lkMatrix);
  arma::rowvec prob_g = Rcpp::as<arma::rowvec>(probg);
  arma::colvec rtab = Rcpp::as<arma::colvec>(ritab);
  arma::colvec wts = Rcpp::as<arma::colvec>(weights);
  arma::colvec awts =  wts * alfa;
  int ng = lkm.n_rows;
  for( int i = 0; i<nrp ; i++){
    arma::mat pggd = arma::diagmat(prob_g) * lkm;
    pggd = pggd * diagmat(1 / sum(pggd,0));
    arma::mat kd = pggd * rtab;
    prob_g = trans(clamp(awts - 1 + kd,0,kd.max()));
    prob_g = prob_g / accu(prob_g);
  }
  return wrap(prob_g);
}
