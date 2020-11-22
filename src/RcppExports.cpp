// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// dmvnrm
arma::vec dmvnrm(arma::mat const& x, arma::rowvec const& mean, arma::mat const& sigma, bool const logd);
RcppExport SEXP _looperr_dmvnrm(SEXP xSEXP, SEXP meanSEXP, SEXP sigmaSEXP, SEXP logdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat const& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::rowvec const& >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< bool const >::type logd(logdSEXP);
    rcpp_result_gen = Rcpp::wrap(dmvnrm(x, mean, sigma, logd));
    return rcpp_result_gen;
END_RCPP
}
// hatdiag
double hatdiag(arma::mat& Q, arma::vec& w, int& i);
RcppExport SEXP _looperr_hatdiag(SEXP QSEXP, SEXP wSEXP, SEXP iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type Q(QSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type w(wSEXP);
    Rcpp::traits::input_parameter< int& >::type i(iSEXP);
    rcpp_result_gen = Rcpp::wrap(hatdiag(Q, w, i));
    return rcpp_result_gen;
END_RCPP
}
// loclin_gauss
Rcpp::List loclin_gauss(arma::mat& X, arma::mat& H, arma::vec& y);
RcppExport SEXP _looperr_loclin_gauss(SEXP XSEXP, SEXP HSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type H(HSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(loclin_gauss(X, H, y));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_looperr_dmvnrm", (DL_FUNC) &_looperr_dmvnrm, 4},
    {"_looperr_hatdiag", (DL_FUNC) &_looperr_hatdiag, 3},
    {"_looperr_loclin_gauss", (DL_FUNC) &_looperr_loclin_gauss, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_looperr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
