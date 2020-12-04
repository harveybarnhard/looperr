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

arma::vec gausskern(arma::mat const &x,
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

// ------------------------------------------------
// Univariate Epanechnikov kernel. Note that there
// are much faster methods, but this more straight
// forward one does the job
// ------------------------------------------------
arma::vec epankern(arma::vec const &x,
                   arma::rowvec const &mean,
                   arma::mat const &h) {
  arma::uword const n = x.n_rows;
  arma::vec out(n, arma::fill::zeros);
  double hsca = arma::as_scalar(h);
  arma::vec u = arma::abs((x.each_row() - mean))/hsca;
  for(arma::uword i = 0; i < n; i++){
    if(arma::as_scalar(u(i))<1){
      out(i) = 3*(1-pow(arma::as_scalar(u(i)), 2))/4;
    }
  }
  return out;
}

// -------------------------------------------------
// Function that finds diagonal of hat matrix given
// the Q of the QR factorization and a vector of
// weights
// -------------------------------------------------
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
//' using Gaussian or Epanechnikov kernel.
//'
//' @param X an nxk data matrix
//' @param H a kxk positive definite bandwidth matrix
//' @param y The nx1 output vector
//' @param Xeval an mxp matrix at which to predict using local linear regression
//' @param sameX binary; Are the evaluation points the same as X? One for yes.
//' @param kernel integer; 1 for Gaussian, 2 for Epanechnikov
// [[Rcpp::export]]
Rcpp::List loclin(arma::mat const &X,
                        arma::mat const &H,
                        arma::vec const &y,
                        arma::mat const &Xeval,
                        int const &sameX,
                        int const &kernel) {
  int n = X.n_rows, k = X.n_cols, neval = Xeval.n_rows;
  arma::uvec colind = arma::regspace<arma::uvec>(1,1,k - 1);
  arma::vec pred_vals(neval, arma::fill::zeros);
  arma::vec hat(n, arma::fill::zeros);
  arma::vec pred_err(n, arma::fill::zeros);
  arma::vec w(n, arma::fill::ones);

  // Perform Cholesky decomposition of H
  arma::mat cholH = H;
  if((k > 2) && (kernel==1)) {
    cholH = arma::chol(H, "lower");
  }
  for(int i=0; i < neval; i++){
    arma::rowvec x0 = Xeval.row(i);
    // Determine weights using a multivariate Gaussian Kernel
    if(kernel==1){
      w = gausskern(X.cols(colind), x0(colind), cholH);
    }
    if(kernel==2){
      w = epankern(X.cols(colind), x0(colind), cholH);
    }
    // Solve OLS using economical QR decomposition, scaling X and y by weights
    arma::mat Q, R;
    arma::qr_econ(Q, R, X.each_col() % sqrt(w));
    arma::vec beta = solve(R, (Q.t() * (y%sqrt(w))));
    pred_vals(i) = arma::as_scalar(x0*beta);
    // Find diagonal of hat matrix
    if(sameX==1){
      hat(i) = hatdiag(Q, w, i);
    }
  }
  if(sameX==1){
    pred_err = (y - pred_vals) / (1 - hat);
  }
  List listout = List::create(Named("fitted.values") = pred_vals,
                              Named("loo_pred_err")  = pred_err,
                              Named("cvscore")       = as_scalar(sum(pow(pred_err,2))));
  return listout;
}

// -----------------------------------------
// Function that performs local linear
// regression for a fixed bandwidth using
// Gaussian kernel
// -----------------------------------------

//' Function that performs local linear regression
//' using Gaussian or Epanechnikov kernel.
//'
//' @param X an nxk data matrix
//' @param H a kxk positive definite bandwidth matrix
//' @param y The nx1 output vector
//' @param Xeval an mxp matrix at which to predict using local linear regression
//' @param sameX binary; Are the evaluation points the same as X? One for yes.
//' @param kernel integer; 1 for Gaussian, 2 for Epanechnikov
//' @param nthr integer; number of threads to use for parallel processing
// [[Rcpp::export]]
Rcpp::List loclin_by(arma::mat const &X,
                  arma::mat const &H,
                  arma::vec const &y,
                  arma::mat const &Xeval,
                  int const &sameX,
                  int const &kernel,
                  int const nthr = 1) {
  int n = X.n_rows, k = X.n_cols, neval = Xeval.n_rows;
  arma::uvec colind = arma::regspace<arma::uvec>(1,1,k - 1);
  arma::vec pred_vals(neval, arma::fill::zeros);
  arma::vec hat(n, arma::fill::zeros);
  arma::vec pred_err(n, arma::fill::zeros);
  arma::vec w(n, arma::fill::ones);

  // Perform Cholesky decomposition of H
  arma::mat cholH = H;
  if((k > 2) && (kernel==1)) {
    cholH = arma::chol(H, "lower");
  }
  for(int i=0; i < neval; i++){
    arma::rowvec x0 = Xeval.row(i);
    // Determine weights using a multivariate Gaussian Kernel
    if(kernel==1){
      w = gausskern(X.cols(colind), x0(colind), cholH);
    }
    if(kernel==2){
      w = epankern(X.cols(colind), x0(colind), cholH);
    }
    // Solve OLS using economical QR decomposition, scaling X and y by weights
    arma::mat Q, R;
    arma::qr_econ(Q, R, X.each_col() % sqrt(w));
    arma::vec beta = solve(R, (Q.t() * (y%sqrt(w))));
    pred_vals(i) = arma::as_scalar(x0*beta);
    // Find diagonal of hat matrix
    if(sameX==1){
      hat(i) = hatdiag(Q, w, i);
    }
  }
  if(sameX==1){
    pred_err = (y - pred_vals) / (1 - hat);
  }
  List listout = List::create(Named("fitted.values") = pred_vals,
                              Named("loo_pred_err")  = pred_err,
                              Named("cvscore")       = as_scalar(sum(pow(pred_err,2))));
  return listout;
}
