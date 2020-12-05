// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <omp.h>
using namespace Rcpp;

// -----------------------------------------
// Function that performs linear regression
// and saves the diagonals of the hat matrix
// -----------------------------------------

//' Function that performs linear regression
//' and saves the diagonal of the hat matrix
//'
//' @param X an nxk numeric data matrix
//' @param y The nx1 numeric output vector
//' @param w an nx1 numeric vector of weights
//' @export
// [[Rcpp::export]]

Rcpp::List fastols(arma::mat const &X, arma::vec const &y, arma::vec const &w) {
  int ncols = X.n_cols;
  // Solve OLS using fast QR decomposition
  arma::vec ws = sqrt(w);
  arma::mat Q, R;
  arma::qr_econ(Q, R, X.each_col() % ws);
  arma::vec Qy(ncols, arma::fill::none);
  for(int k = 0; k < ncols; k++){
    Qy(k) = dot(Q.col(k), y % ws);
  }
  arma::vec beta = solve(R, Qy);
  // Find diagonal of hat matrix given the Q of the QR factorization above
  arma::mat Q1 = Q.each_col() / ws;
  arma::mat Q2 = Q.each_col() % ws;
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

// -----------------------------------------
// Same as fastols above but over groups
// -----------------------------------------

//' Function that performs linear regression
//' and saves the diagonal of the hat matrix
//'
//' @param X an nxk numeric data matrix
//' @param y The nx1 numeric output vector
//' @param w an nx1 numeric vector of weights
//' @param g an nx1 sorted integer vector of groups
//' @param nthr integer; number of threads to use for parallel processing
// [[Rcpp::export]]
Rcpp::List fastols_by2(arma::mat const &X,
                      arma::vec const &y,
                      arma::vec const &w,
                      IntegerVector const &g,
                      int const nthr = 1,
                      int const compute_se = 1,
                      int const compute_hat = 0) {
  int nrows = X.n_rows, ncols = X.n_cols;
  arma::vec beta(nrows, arma::fill::zeros),
  hat(nrows, arma::fill::zeros),
  pred_err(nrows, arma::fill::zeros),
  groups(nrows, arma::fill::zeros),
  start(nrows, arma::fill::zeros),
  end(nrows, arma::fill::zeros);
  omp_set_num_threads(nthr);
  // Weight
  arma::vec ws = sqrt(w);
  arma::mat Xw = X.each_col() % ws;
  arma::vec yw = y % ws;

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
  // Initialize matrix to store coefficients and loop over groups
  arma::mat betamat(ncols, numgrps, arma::fill::none);
  arma::mat varmat(ncols, numgrps, arma::fill::zeros);
#pragma omp parallel for
  for(int j=0; j < numgrps; j++){
    // Solve OLS using fast QR decomposition
    arma::mat Q, R;
    arma::qr_econ(Q, R, Xw.rows(start(j), end(j)));
    arma::vec Qy(ncols, arma::fill::none);
    for(int k = 0; k < ncols; k++){
      Qy(k) = dot(Q.col(k), yw.subvec(start(j), end(j)));
    }
    betamat.col(j) = solve(R, Qy);
    // Find residuals noting that Xbeta = QQ'y
    yw.subvec(start(j), end(j)) = yw.subvec(start(j), end(j)) - Q*Qy;
    // Find standard OLS variance
    if(compute_se == 1){
      double sig2 = arma::as_scalar(arma::sum(pow(yw.subvec(start(j), end(j)), 2))/(Q.n_rows - ncols));
      R = inv(R);
      varmat.col(j) = sig2*arma::sum(R % R, 1);
    }
    // Find diagonal of hat matrix given the Q of the QR factorization above
    if(compute_hat==1){
      hat.subvec(start(j), end(j)) = sum( (Q.each_col() / ws.subvec(start(j), end(j))) %
        (Q.each_col() % ws.subvec(start(j), end(j))), 1);
    }
  }
  // Calculate leave-one-out prediction error
  List listout = List::create(Named("beta")          = betamat,
                              Named("hatdiag")       = hat,
                              Named("residuals")     = yw,
                              Named("variance")      = varmat);
  return listout;
}
