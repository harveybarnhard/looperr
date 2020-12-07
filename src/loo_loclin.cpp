// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <omp.h>
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

// ------------------------------------------------
// Uniform kernel
// ------------------------------------------------
arma::vec unifkern(arma::vec const &x,
                   arma::rowvec const &mean,
                   arma::mat const &h) {
  arma::uword const n = x.n_rows;
  arma::vec out(n, arma::fill::zeros);
  double hsca = arma::as_scalar(h);
  arma::vec u = arma::abs((x.each_row() - mean))/hsca;
  for(arma::uword i = 0; i < n; i++){
    if(arma::as_scalar(u(i))<1){
      out(i) = 1/(2*hsca);
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
                        int const &kernel,
                        int const nthr = 1) {
  int nrows = X.n_rows, ncols = X.n_cols, neval = Xeval.n_rows;
  arma::uvec colind = arma::regspace<arma::uvec>(1,1,ncols - 1);
  arma::vec pred_vals(neval, arma::fill::zeros);
  arma::vec ws(nrows, arma::fill::ones);
  omp_set_num_threads(nthr);
  // Perform Cholesky decomposition of H
  arma::mat cholH = H;
  if((ncols > 2) && (kernel==1)) {
    cholH = arma::chol(H, "lower");
  }
#pragma omp parallel for
  for(int i=0; i < neval; i++){
    arma::rowvec x0 = Xeval.row(i);
    // Determine weights using a multivariate Gaussian Kernel
    if(kernel==1){
      ws = sqrt(gausskern(X.cols(colind), x0(colind), cholH));
    }
    else if(kernel==2){
      ws = sqrt(epankern(X.cols(colind), x0(colind), cholH));
    }
    else if(kernel==3){
      ws = sqrt(unifkern(X.cols(colind), x0(colind), cholH));
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
                        int const &kernel,
                        int const nthr = 1) {
  int nrows = X.n_rows, ncols = X.n_cols;
  arma::uvec colind = arma::regspace<arma::uvec>(1,1,ncols - 1);
  arma::vec pred_vals(nrows, arma::fill::zeros);
  arma::vec hat(nrows, arma::fill::zeros);
  arma::vec pred_err(nrows, arma::fill::zeros);
  arma::vec ws(nrows, arma::fill::ones);
  omp_set_num_threads(nthr);

  // Perform Cholesky decomposition of H
  arma::mat cholH = H;
  if((ncols > 2) && (kernel==1)) {
    cholH = arma::chol(H, "lower");
  }
#pragma omp parallel for
  for(int i=0; i < nrows; i++){
    arma::rowvec x0 = X.row(i);
    // Determine weights using a multivariate Gaussian Kernel
    if(kernel==1){
      ws = sqrt(gausskern(X.cols(colind), x0(colind), cholH));
    }
    else if(kernel==2){
      ws = sqrt(epankern(X.cols(colind), x0(colind), cholH));
    }
    else if(kernel==3){
      ws = sqrt(unifkern(X.cols(colind), x0(colind), cholH));
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
//' @param y The nx1 output vector
//' @param H a kxk positive definite bandwidth matrix
//' @param kernel integer; 1 for Gaussian, 2 for Epanechnikov
// [[Rcpp::export]]
Rcpp::List loclin_sameX2(arma::mat const &X,
                        arma::vec const &y,
                        double const &H,
                        int const &kernel,
                        int const nthr = 1) {
  int nrows = X.n_rows, ncols = X.n_cols;
  arma::uvec colind = arma::regspace<arma::uvec>(1,1,ncols - 1);
  arma::vec pred_vals(nrows, arma::fill::zeros);
  arma::vec hat(nrows, arma::fill::zeros);
  arma::vec pred_err(nrows, arma::fill::zeros);
  arma::vec ws(nrows, arma::fill::ones);
  omp_set_num_threads(nthr);

  // Initialize first and last element for compact kernels
  int jlast=0, jfirst=0;

  for(int i=0; i < nrows; i++){
    arma::rowvec x0 = X.row(i);
    // Determine weights using a multivariate Gaussian Kernel
    if(kernel==2){
      // First find the new lower end
      for(int j=jfirst; j < nrows; j++){
        double test = arma::as_scalar(x0(colind) -  X.cols(colind));
        if(test > H){
          jfirst = j - 1;
          break;
        }
      }
      // Find the new upper end
      for(int j=jlast; j < nrows; j++){
        double test = arma::as_scalar(X.cols(colind) - x0(colind));
        if(test > H){
          jlast = j - 1;
          break;
        }
      }
      arma::vec u = X.submat(jfirst, colind, jlast, colind);
      arma::vec u = arma::abs(arma::colvec X.submat(jfirst, colind, jlast, colind) - x0(colind))/H;
      arma::vec ws = sqrt(3*(1-pow(u, 2))/4);
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
*/

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
Rcpp::List loclin_sameX2(arma::mat const &X,
                         arma::vec const &y,
                         double const &H,
                         int const &kernel,
                         int const nthr = 1) {
  int nrows = X.n_rows, ncols = X.n_cols;
  arma::vec pred_vals(nrows, arma::fill::zeros);
  arma::vec hat(nrows, arma::fill::zeros);
  // Initialize first and last element for compact kernels
  int jfirst=0, jlast=0;

  for(int i=0; i < nrows; i++){
    Rcout << "Iteration " << i + 1 << std::endl;
    // First find the new lower bound of kernel support
    while(jfirst < nrows - 1){
      double test = arma::as_scalar(X(i,1) - X(jfirst,1));
      if(test < H){
        break;
      }
      jfirst++;
    }
    // Find the new upper bound of kernel support
    while(jlast < nrows - 1){
      double test = arma::as_scalar(X(jlast, 1) - X(i,1));
      if(test > H){
        jlast--;
        break;
      }
      jlast++;
    }
    if(jfirst==jlast){
      pred_vals(i) = arma::datum::nan;
      hat(i) = arma::datum::nan;
      continue;
    }
    arma::vec u = arma::abs(X.submat(jfirst, 1, jlast, 1) - X(i,1));
    arma::vec ws = sqrt(3*(1-pow(u/H, 2))/4);
    arma::mat Xw = X.rows(jfirst, jlast);
    arma::mat Q, R;
    arma::qr_econ(Q, R, Xw.each_col() % ws);
    Rcout << "Unscaled prediction " << 1-pow(u/H, 2) << std::endl;
    arma::vec Qy(ncols, arma::fill::none);
    arma::vec yw = y.subvec(jfirst, jlast) % ws;
    for(int k = 0; k < ncols; k++){
      Qy(k) = dot(Q.col(k), yw);
    }
    // Find fitted values and the ith diagonal of hat matrix. Note that
    // beta does not need to be computed for same X because
    // yhat = Xbetahat = QQ'y
    int selfind = i - jfirst;
    pred_vals(i) = arma::as_scalar((Q.row(selfind)*Qy)/ws(selfind));
    hat(i) = dot(Q.row(selfind), Q.row(selfind));
  }
  List listout = List::create(Named("hat") = hat,
                              Named("predvals") = pred_vals);
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
//' @param kernel integer; 1 for Gaussian, 2 for Epanechnikov, 3 for uniform
//' @param nthr integer; number of threads to use with OpenMP
// [[Rcpp::export]]
Rcpp::List loclin_sameX_by(arma::mat const &X,
                        arma::vec const &y,
                        IntegerVector const &g,
                        arma::mat const &H,
                        int const &kernel,
                        int const nthr = 1) {
  int nrows = X.n_rows, ncols = X.n_cols;
  arma::uvec colind = arma::regspace<arma::uvec>(1,1,ncols - 1);
  arma::vec pred_vals(nrows, arma::fill::zeros);
  arma::vec hat(nrows, arma::fill::zeros);
  arma::vec pred_err(nrows, arma::fill::zeros);
  arma::vec ws(nrows, arma::fill::ones);
  omp_set_num_threads(nthr);

  // Perform Cholesky decomposition of H
  arma::mat cholH = H;
  if((ncols > 2) && (kernel==1)) {
    cholH = arma::chol(H, "lower");
  }
#pragma omp parallel for
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

