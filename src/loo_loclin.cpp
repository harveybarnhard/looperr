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

// -----------------------------------------
// Function that performs local linear
// regression for a fixed bandwidth using
// Gaussian kernel
// -----------------------------------------

//' Function that performs local linear regression
//' using Gaussian or Epanechnikov kernel.
//' @param X an nxk data matrix
//' @param y The nx1 output vector
//' @param H a kxk positive definite bandwidth matrix
//' @param Xeval an mxp matrix at which to predict using local linear regression
//' @param kernel integer; 1 for Gaussian, 2 for Epanechnikov
// [[Rcpp::export]]
Rcpp::List loclin_diffX(arma::mat const &X,
                        arma::vec const &y,
                        arma::mat const &H,
                        arma::mat const &Xeval,
                        int const &kernel) {
  int nrows = X.n_rows, ncols = X.n_cols, neval = Xeval.n_rows;
  arma::uvec colind = arma::regspace<arma::uvec>(1,1,ncols - 1);
  arma::vec pred_vals(neval, arma::fill::zeros);
  arma::vec ws(nrows, arma::fill::ones);

  // Perform Cholesky decomposition of H
  arma::mat cholH = H;
  if((ncols > 2) && (kernel==1)) {
    cholH = arma::chol(H, "lower");
  }
  for(int i=0; i < neval; i++){
    arma::rowvec x0 = Xeval.row(i);
    // Determine weights using a multivariate Gaussian Kernel
    if(kernel==1){
      ws = sqrt(gausskern(X.cols(colind), x0(colind), cholH));
    }
    if(kernel==2){
      ws = sqrt(epankern(X.cols(colind), x0(colind), cholH));
    }
    // Solve OLS using economical QR decomposition, scaling X and y by weights
    arma::mat Q, R;
    arma::qr_econ(Q, R, X.each_col() % ws);
    arma::vec Qy(ncols, arma::fill::none);
    arma::vec yw = y % ws;
    for(int k = 0; k < ncols; k++){
      Qy(k) = dot(Q.col(k), yw);
    }
    arma::vec beta = solve(R, Qy);
    pred_vals(i) = arma::as_scalar(x0*beta);
  }
  List listout = List::create(Named("fitted.values") = pred_vals);
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
//' @param y The nx1 output vector
//' @param H a kxk positive definite bandwidth matrix
//' @param kernel integer; 1 for Gaussian, 2 for Epanechnikov
// [[Rcpp::export]]
Rcpp::List loclin_sameX(arma::mat const &X,
                        arma::vec const &y,
                        arma::mat const &H,
                        int const &kernel) {
  int nrows = X.n_rows, ncols = X.n_cols;
  arma::uvec colind = arma::regspace<arma::uvec>(1,1,ncols - 1);
  arma::vec pred_vals(nrows, arma::fill::zeros);
  arma::vec hat(nrows, arma::fill::zeros);
  arma::vec pred_err(nrows, arma::fill::zeros);
  arma::vec ws(nrows, arma::fill::ones);

  // Perform Cholesky decomposition of H
  arma::mat cholH = H;
  if((ncols > 2) && (kernel==1)) {
    cholH = arma::chol(H, "lower");
  }
  for(int i=0; i < nrows; i++){
    arma::rowvec x0 = X.row(i);
    // Determine weights using a multivariate Gaussian Kernel
    if(kernel==1){
      ws = sqrt(gausskern(X.cols(colind), x0(colind), cholH));
    }
    if(kernel==2){
      ws = sqrt(epankern(X.cols(colind), x0(colind), cholH));
    }
    // Solve OLS using economical QR decomposition, scaling X and y by weights
    arma::mat Q, R;
    arma::qr_econ(Q, R, X.each_col() % ws);
    arma::vec Qy(ncols, arma::fill::none);
    arma::vec yw = y % ws;
    for(int k = 0; k < ncols; k++){
      Qy(k) = dot(Q.col(k), yw);
    }
    // Find fitted values and the ith diagonal of hat matrix. Note that
    // beta does not need to be computed for same X because
    // yhat = Xbetahat = QQ'y
    pred_vals(i) = arma::as_scalar((Q.row(i)*Qy)/ws(i));
    hat(i) = dot(Q.row(i), Q.row(i));
  }
  pred_err = (y - pred_vals)/(1-hat);
  List listout = List::create(Named("fitted.values") = pred_vals,
                              Named("loo_pred_err")  = pred_err,
                              Named("cvscore")       = arma::dot(pred_err, pred_err));
  return listout;
}
/*
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
                  IntegerVector const &g,
                  arma::mat const &Xeval,
                  int const &sameX,
                  int const &kernel,
                  int const nthr = 1) {
  int nrows = X.n_rows, ncols = X.n_cols, neval = Xeval.n_rows;
  arma::uvec colind = arma::regspace<arma::uvec>(1,1,ncols - 1);
  arma::vec pred_vals(neval, arma::fill::zeros),
  hat(nrows, arma::fill::zeros),
  pred_err(nrows, arma::fill::zeros),
  ws(nrows, arma::fill::ones),
  start(nrows, arma::fill::zeros),
  end(nrows, arma::fill::zeros);
  omp_set_num_threads(nthr);

  // Perform Cholesky decomposition of H
  arma::mat cholH = H;
  if((ncols > 2) && (kernel==1)) {
    cholH = arma::chol(H, "lower");
  }

  // Find start and endpoints of group, assuming g is already sorted
  // Find start and endpoints of group, assuming g is already sorted
  int cur = g(0);
  int numgrps = 1;
  for(int i=0; i < nrows; i++){
    if(g(i)!=cur){
      end(numgrps - 1) = i - 1;
      numgrps += 1;
      start(numgrps - 1) = i;
      cur = g(i);
    }
  }
  end(numgrps - 1) = nrows - 1;
  start = start.head(numgrps);
  end   = end.head(numgrps);
#pragma omp parallel for
  for(int i=0; i < neval; i++){
    arma::rowvec x0 = Xeval.row(i);
    // Determine weights using a multivariate Gaussian Kernel
    if(kernel==1){
      ws = sqrt(gausskern(X.cols(colind), x0(colind), cholH));
    }
    if(kernel==2){
      ws = sqrt(epankern(X.cols(colind), x0(colind), cholH));
    }
    // Solve OLS using economical QR decomposition, scaling X and y by weights
    arma::mat Q, R;
    arma::qr_econ(Q, R, X.each_col() % ws);
    arma::vec Qy(ncols, arma::fill::none);
    arma::vec yw = y % ws;
    for(int k = 0; k < nrows; k++){
      Qy(k) = dot(Q.col(k), yw);
    }
    arma::vec beta = solve(R, Qy);
    // TODO Optimize prediction to only predict ONE value
    pred_vals(i) = arma::as_scalar(Q*Qy/ws);
    // Find diagonal of hat matrix
    if(sameX==1){
      hat(i) = hatdiag(Q, pow(ws,2), i);
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
*/
