// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// fastols
Rcpp::List fastols(arma::mat const& X, arma::vec const& y, arma::vec const& w);
RcppExport SEXP _looperr_fastols(SEXP XSEXP, SEXP ySEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat const& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec const& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec const& >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(fastols(X, y, w));
    return rcpp_result_gen;
END_RCPP
}
// fastols_by
Rcpp::List fastols_by(arma::mat const& X, arma::vec const& y, arma::vec const& w, IntegerVector const& g, int const nthr);
RcppExport SEXP _looperr_fastols_by(SEXP XSEXP, SEXP ySEXP, SEXP wSEXP, SEXP gSEXP, SEXP nthrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat const& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec const& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec const& >::type w(wSEXP);
    Rcpp::traits::input_parameter< IntegerVector const& >::type g(gSEXP);
    Rcpp::traits::input_parameter< int const >::type nthr(nthrSEXP);
    rcpp_result_gen = Rcpp::wrap(fastols_by(X, y, w, g, nthr));
    return rcpp_result_gen;
END_RCPP
}
// loclin
Rcpp::List loclin(arma::mat const& X, arma::mat const& H, arma::vec const& y, arma::mat const& Xeval, int const& sameX, int const& kernel);
RcppExport SEXP _looperr_loclin(SEXP XSEXP, SEXP HSEXP, SEXP ySEXP, SEXP XevalSEXP, SEXP sameXSEXP, SEXP kernelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat const& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type H(HSEXP);
    Rcpp::traits::input_parameter< arma::vec const& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type Xeval(XevalSEXP);
    Rcpp::traits::input_parameter< int const& >::type sameX(sameXSEXP);
    Rcpp::traits::input_parameter< int const& >::type kernel(kernelSEXP);
    rcpp_result_gen = Rcpp::wrap(loclin(X, H, y, Xeval, sameX, kernel));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_looperr_fastols", (DL_FUNC) &_looperr_fastols, 3},
    {"_looperr_fastols_by", (DL_FUNC) &_looperr_fastols_by, 5},
    {"_looperr_loclin", (DL_FUNC) &_looperr_loclin, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_looperr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
