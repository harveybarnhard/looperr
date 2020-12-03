// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;

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
//' @param X: an nxk data matrix
//' @param y: The nx1 output vector
//' @param w: an nx1 vector of weights
// [[Rcpp::export]]
Rcpp::List fastols(arma::mat const &X,
                   arma::vec const &y,
                   arma::vec const &w,
                   IntegerVector const &g) {
  int nrows = X.n_rows;
  arma::vec beta(nrows, arma::fill::zeros),
            hat(nrows, arma::fill::zeros),
            pred_err(nrows, arma::fill::zeros),
            fittedy(nrows, arma::fill::zeros),
            groups(nrows, arma::fill::zeros),
            start(nrows, arma::fill::zeros),
            end(nrows, arma::fill::zeros);
  // Find start and endpoints of groupings, assuming g is already sorted
  int cur = g(0);
  int numgrps = 0;
  for(int i=0; i < nrows; i++){
    if(g(i)!=cur){
      end(numgrps) = i - 1;
      numgrps += 1;
      start(numgrps) = i;
      cur = g(i);
    }
  }
  end(numgrps) = nrows;
  start = start.head(numgrps);
  end   = end.head(numgrps);

  for(int j=0; j < numgrps; j++){
    // Solve OLS using fast QR decomposition
    arma::mat Q, R;
    arma::qr_econ(Q, R, X.row(start(j)).each_col() % sqrt(w));
    arma::vec beta = solve(R, (Q.t()*(y % sqrt(w))));
    // Find diagonal of hat matrix given the Q of the QR factorization above
    arma::mat Q1 = Q.each_col() / sqrt(w);
    arma::mat Q2 = Q.each_col() % sqrt(w);
    arma::vec hat = sum(Q1%Q2,1);
    // Find LOO prediction error
    arma::vec fittedy = X*beta;
    arma::vec pred_err = (y - fittedy) / (1 - hat);
  }
  List listout = List::create(Named("beta")          = beta,
                              Named("hatdiag")       = hat,
                              Named("loo_pred_err")  = pred_err,
                              Named("fitted.values") = fittedy);
  return listout;
}
