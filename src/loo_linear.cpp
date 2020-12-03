// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
# include <omp.h>
using namespace Rcpp;
using namespace arma;

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
Rcpp::List fastols_by(arma::mat const &X,
                   arma::vec const &y,
                   arma::vec const &w,
                   IntegerVector const &g,
                   int const nthr = 1) {
  int nrows = X.n_rows;
  arma::vec beta(nrows, arma::fill::zeros),
            hat(nrows, arma::fill::zeros),
            pred_err(nrows, arma::fill::zeros),
            fittedy(nrows, arma::fill::zeros),
            groups(nrows, arma::fill::zeros),
            start(nrows, arma::fill::zeros),
            end(nrows, arma::fill::zeros);
  omp_set_num_threads(nthr);
  // Find start and endpoints of groupings, assuming g is already sorted
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
  arma::mat betamat(X.n_cols, numgrps, arma::fill::zeros);
#pragma omp parallel for
  for(int j=0; j < numgrps; j++){
    // Solve OLS using fast QR decomposition
    arma::mat Xj = X.rows(start(j), end(j));
    arma::vec yj = y.subvec(start(j), end(j));
    arma::vec wj = w.subvec(start(j), end(j));
    arma::mat Q, R;
    arma::qr_econ(Q, R, Xj.each_col() % sqrt(wj));
    arma::vec beta = solve(R, (Q.t()*(yj % sqrt(wj))));
    betamat.col(j) = beta;
    // Find diagonal of hat matrix given the Q of the QR factorization above
    arma::mat Q1 = Q.each_col() / sqrt(wj);
    arma::mat Q2 = Q.each_col() % sqrt(wj);
    hat.subvec(start(j), end(j)) = sum(Q1%Q2,1);
    // Find fitted y values
    fittedy.subvec(start(j), end(j)) = Xj*beta;
  }
  // Calculate leave-one-out prediction error
  pred_err = (y - fittedy) / (1 - hat);
  List listout = List::create(Named("beta")          = betamat,
                              Named("hatdiag")       = hat,
                              Named("loo_pred_err")  = pred_err,
                              Named("fitted.values") = fittedy);
  return listout;
}
