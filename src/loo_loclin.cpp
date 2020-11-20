#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// -----------------------------
// Multivariate density function: copied from
// https://gallery.rcpp.org/articles/dmvnorm_arma/
// -----------------------------
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
                           arma::mat const &sigma,
                           bool const logd = false) {
  using arma::uword;
  uword const n = x.n_rows,
    xdim = x.n_cols;
  arma::vec out(n);
  arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
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


// ---------------------------------------
// Function that performs OLS regression template using
// https://dannyjameswilliams.co.uk/portfolios/sc2/rcpp/
// ---------------------------------------

// [[Rcpp::export(name="fastlm")]]
arma::vec fastols(arma::vec& y, arma::mat& X) {
  // Solve OLS using fast QR decomposition
  arma::mat Q, R;
  arma::qr_econ(Q, R, X);
  arma::vec beta = solve(R, (Q.t() * y));
  return beta;
}

// ------------------------
// Function that performs one
// ------------------------
// [[Rcpp::export(name="testh")]]
arma::vec loclin_row(arma::rowvec& x0,
                     arma::mat& X,
                     arma::mat& H) {
  // Determine weights using a multivariate Gaussian Kernel
  arma::vec h = dmvnrm(X, x0.t(), H);
  /*
  arma::vec beta = fastols(y % sqrt(h), X.each_col() % sqrt(h));
  return X0 * beta;
  */
}

/*** R

*/
