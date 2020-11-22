// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

// ------------------------------------------------
// Multivariate density function: modified from
// https://gallery.rcpp.org/articles/dmvnorm_arma/
// ------------------------------------------------
static double const log2pi = std::log(2.0 * M_PI);
/* C++ version of the dtrmv BLAS function */
void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat){
  arma::uword const n = trimat.n_cols;

  for(unsigned j = n; j-- > 0;){
    double tmp(0.);
    for(unsigned i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x[i];
    x[j] = tmp;
  }
}

// [[Rcpp::export]]
arma::vec dmvnrm(arma::mat const &x,
                           arma::rowvec const &mean,
                           arma::mat const &chol_sigma,
                           bool const logd = false) {
  using arma::uword;
  uword const n = x.n_rows,
    xdim = x.n_cols;
  arma::vec out(n);
  arma::mat const rooti = arma::inv(trimatu(chol_sigma));
  double const rootisum = arma::sum(log(rooti.diag())),
    constants = -(double)xdim/2.0 * log2pi,
    other_terms = rootisum + constants;

  arma::rowvec z;
  for (uword i = 0; i < n; i++) {
    z = (x.row(i) - mean);
    inplace_tri_mat_mult(z, rooti);
    out(i) = other_terms - 0.5 * arma::dot(z, z);
  }

  if (logd)
    return out;
  return exp(out);
}


// ------------------------------------------------------
// Function that performs OLS regression template using
// https://dannyjameswilliams.co.uk/portfolios/sc2/rcpp/
// ------------------------------------------------------
arma::vec fastols(arma::vec const &y, arma::mat const &X) {
  // Solve OLS using fast QR decomposition
  arma::mat Q, R;
  arma::qr_econ(Q, R, X);
  arma::vec beta = solve(R, (Q.t() * y));
  return beta;
}

// -------------------------------------------------
// Function that finds diagonal of hat matrix given
// the Q of the QR factorization and a vector of
// weights
// https://stackoverflow.com/questions/20562177/get-hat-matrix-from-qr-decomposition-for-weighted-least-square-regression
// -------------------------------------------------
// [[Rcpp::export]]
double hatdiag(arma::mat const &Q, arma::vec const &w, int const &i) {
  // Rescale Q by weights
  arma::mat Q1 = Q.each_col() / sqrt(w);
  arma::mat Q2 = Q.each_col() % sqrt(w);
  arma::vec out = sum(Q1%Q2,1);
  return out(i);
}

// -----------------------------------------
// Function that performs local linear
// regression for a fixed bandwidth using
// Gaussian kernel
// -----------------------------------------

//' Function that performs local linear regression
//' using a Gaussian Kernel.
//'
//' @param \code{X}: an nxk data matrix
//' @param \code{H}: a kxk positive definite bandwidth matrix
//' @param \code{y}: The nx1 output vector
//' @export
// [[Rcpp::export]]
Rcpp::List loclin_gauss(arma::mat& X, arma::mat& H, arma::vec& y) {
  int n = X.n_rows, k = X.n_cols;
  arma::uvec colind = arma::regspace<arma::uvec>(1,1,k - 1);
  arma::vec pred_vals(n, arma::fill::zeros);
  arma::vec hat(n, arma::fill::zeros);
  arma::rowvec x0(k, arma::fill::zeros);

  // Perform Cholesky decomposition of H
  arma::mat cholH = H;
  if(k > 2) {
    cholH = arma::chol(H);
  }
  for(int i=0; i < n; i++){
    x0 = X.row(i);
    // Determine weights using a multivariate Gaussian Kernel
    arma::vec w = dmvnrm(X.cols(colind), x0(colind), cholH);
    // Scale X and Y by weights
    // X.each_col() %= sqrt(w);
    // Solve OLS using economical QR decomposition
    arma::mat Q, R;
    arma::qr_econ(Q, R, X.each_col() % sqrt(w));
    arma::vec beta = solve(R, (Q.t() * (y%sqrt(w))));
    pred_vals(i) = arma::as_scalar(x0*beta);
    // Find diagonal of hat matrix
    hat(i) = hatdiag(Q, w, i);
  }
  List listout = List::create(Named("pred_vals") = pred_vals,
                              Named("hatdiag") = hat);
  return listout;
}
//' Function that returns LOOCV score using Gaussian kernel
//'
//' @param \code{X}: an nxk data matrix
//' @param \code{H}: a kxk positive definite bandwidth matrix
//' @param \code{y}: The nx1 output vector
//' @export
// [[Rcpp::export]]
double loocv(arma::mat& X, arma::mat& H, arma::vec& y) {
  Rcpp::List L = loclin_gauss(X, H, y);
  arma::vec pred_vals = L["pred_vals"];
  arma::vec hatdiag = L["hatdiag"];
  double cvscore = sum(pow((y - pred_vals) / (1 - hatdiag),2));
  return cvscore;
}
