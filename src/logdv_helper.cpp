#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export(rng = FALSE)]]
arma::vec diag_cpp(arma::mat & y, arma::mat & S){
  // compute diagonal
  arma::mat M = y * S * y.t();

  return M.diag();

  // update candidates
  // Rcpp::Rcout << "m = " << m << "\n";
  // Rcpp::Rcout << "one-by-one = " << arma::vec(1).fill(index) << "\n";
  //return arma::conv_to<std::vector<double>>::from(K);

}
