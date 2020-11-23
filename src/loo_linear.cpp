// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

// -----------------------------------------
// Function that performs linear regression
// and saves the diagonals of the hat matrix
// -----------------------------------------

//' Function that performs linear regression
//' and saves the diagonal of the hat matrix
//'
//' @param X: an nxk data matrix
//' @param y: The nx1 output vector
//' @param w: an nx1 vector of weights
//' @export
// [[Rcpp::export]]
Rcpp::List fastols(arma::mat& X, arma::vec& y, arma::vec& w) {
  // Solve OLS using fast QR decomposition
  arma::mat Q, R;
  arma::qr_econ(Q, R, X.each_col() % sqrt(w));
  arma::vec beta = solve(R, (Q.t()*(y % sqrt(w))));
  // Find diagonal of hat matrix given the Q of the QR factorization above
  arma::mat Q1 = Q.each_col() / sqrt(w);
  arma::mat Q2 = Q.each_col() % sqrt(w);
  arma::vec hat = sum(Q1%Q2,1);
  // Find LOO prediction error
  arma::vec fittedy = X*beta;
  arma::vec pred_err = (y - fittedy) / (1 - hat);
  List listout = List::create(Named("beta")          = beta,
                              Named("hatdiag")       = hat,
                              Named("loo_pred_err")  = pred_err,
                              Named("fitted.values") = fittedy);
  return listout;
}
