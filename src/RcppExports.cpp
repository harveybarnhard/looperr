// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// fastols
Rcpp::List fastols(arma::mat const& X, arma::vec const& y, int const compute_se, int const compute_hat);
RcppExport SEXP _looperr_fastols(SEXP XSEXP, SEXP ySEXP, SEXP compute_seSEXP, SEXP compute_hatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat const& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec const& >::type y(ySEXP);
    Rcpp::traits::input_parameter< int const >::type compute_se(compute_seSEXP);
    Rcpp::traits::input_parameter< int const >::type compute_hat(compute_hatSEXP);
    rcpp_result_gen = Rcpp::wrap(fastols(X, y, compute_se, compute_hat));
    return rcpp_result_gen;
END_RCPP
}
// fastolswt
Rcpp::List fastolswt(arma::mat const& X, arma::vec const& y, arma::vec const& w, int const compute_se, int const compute_hat);
RcppExport SEXP _looperr_fastolswt(SEXP XSEXP, SEXP ySEXP, SEXP wSEXP, SEXP compute_seSEXP, SEXP compute_hatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat const& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec const& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec const& >::type w(wSEXP);
    Rcpp::traits::input_parameter< int const >::type compute_se(compute_seSEXP);
    Rcpp::traits::input_parameter< int const >::type compute_hat(compute_hatSEXP);
    rcpp_result_gen = Rcpp::wrap(fastolswt(X, y, w, compute_se, compute_hat));
    return rcpp_result_gen;
END_RCPP
}
// fastols_by
Rcpp::List fastols_by(arma::mat const& X, arma::vec const& y, IntegerVector const& g, int const nthr, int const compute_se, int const compute_hat);
RcppExport SEXP _looperr_fastols_by(SEXP XSEXP, SEXP ySEXP, SEXP gSEXP, SEXP nthrSEXP, SEXP compute_seSEXP, SEXP compute_hatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat const& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec const& >::type y(ySEXP);
    Rcpp::traits::input_parameter< IntegerVector const& >::type g(gSEXP);
    Rcpp::traits::input_parameter< int const >::type nthr(nthrSEXP);
    Rcpp::traits::input_parameter< int const >::type compute_se(compute_seSEXP);
    Rcpp::traits::input_parameter< int const >::type compute_hat(compute_hatSEXP);
    rcpp_result_gen = Rcpp::wrap(fastols_by(X, y, g, nthr, compute_se, compute_hat));
    return rcpp_result_gen;
END_RCPP
}
// fastols_bywt
Rcpp::List fastols_bywt(arma::mat const& X, arma::vec const& y, arma::vec const& w, IntegerVector const& g, int const nthr, int const compute_se, int const compute_hat);
RcppExport SEXP _looperr_fastols_bywt(SEXP XSEXP, SEXP ySEXP, SEXP wSEXP, SEXP gSEXP, SEXP nthrSEXP, SEXP compute_seSEXP, SEXP compute_hatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat const& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec const& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec const& >::type w(wSEXP);
    Rcpp::traits::input_parameter< IntegerVector const& >::type g(gSEXP);
    Rcpp::traits::input_parameter< int const >::type nthr(nthrSEXP);
    Rcpp::traits::input_parameter< int const >::type compute_se(compute_seSEXP);
    Rcpp::traits::input_parameter< int const >::type compute_hat(compute_hatSEXP);
    rcpp_result_gen = Rcpp::wrap(fastols_bywt(X, y, w, g, nthr, compute_se, compute_hat));
    return rcpp_result_gen;
END_RCPP
}
// loclin_diffX
Rcpp::List loclin_diffX(arma::mat const& X, arma::vec const& y, arma::mat const& H, arma::mat const& Xeval, int const& kernel, int const nthr);
RcppExport SEXP _looperr_loclin_diffX(SEXP XSEXP, SEXP ySEXP, SEXP HSEXP, SEXP XevalSEXP, SEXP kernelSEXP, SEXP nthrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat const& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec const& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type H(HSEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type Xeval(XevalSEXP);
    Rcpp::traits::input_parameter< int const& >::type kernel(kernelSEXP);
    Rcpp::traits::input_parameter< int const >::type nthr(nthrSEXP);
    rcpp_result_gen = Rcpp::wrap(loclin_diffX(X, y, H, Xeval, kernel, nthr));
    return rcpp_result_gen;
END_RCPP
}
// loclin_sameX
Rcpp::List loclin_sameX(arma::mat const& X, arma::vec const& y, arma::mat const& H, int const& kernel, int const nthr);
RcppExport SEXP _looperr_loclin_sameX(SEXP XSEXP, SEXP ySEXP, SEXP HSEXP, SEXP kernelSEXP, SEXP nthrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat const& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec const& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type H(HSEXP);
    Rcpp::traits::input_parameter< int const& >::type kernel(kernelSEXP);
    Rcpp::traits::input_parameter< int const >::type nthr(nthrSEXP);
    rcpp_result_gen = Rcpp::wrap(loclin_sameX(X, y, H, kernel, nthr));
    return rcpp_result_gen;
END_RCPP
}
// loclin_sameX_unif
Rcpp::List loclin_sameX_unif(arma::mat const& X, arma::vec const& y, double const& H);
RcppExport SEXP _looperr_loclin_sameX_unif(SEXP XSEXP, SEXP ySEXP, SEXP HSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat const& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec const& >::type y(ySEXP);
    Rcpp::traits::input_parameter< double const& >::type H(HSEXP);
    rcpp_result_gen = Rcpp::wrap(loclin_sameX_unif(X, y, H));
    return rcpp_result_gen;
END_RCPP
}
// loclin_diffX_unif
Rcpp::List loclin_diffX_unif(arma::mat const& X, arma::vec const& y, double const& H, arma::mat const& Xeval);
RcppExport SEXP _looperr_loclin_diffX_unif(SEXP XSEXP, SEXP ySEXP, SEXP HSEXP, SEXP XevalSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat const& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec const& >::type y(ySEXP);
    Rcpp::traits::input_parameter< double const& >::type H(HSEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type Xeval(XevalSEXP);
    rcpp_result_gen = Rcpp::wrap(loclin_diffX_unif(X, y, H, Xeval));
    return rcpp_result_gen;
END_RCPP
}
// loclin_sameX_unif_by
Rcpp::List loclin_sameX_unif_by(arma::mat const& X, arma::vec const& y, IntegerVector const& g, double const& H, int const nthr);
RcppExport SEXP _looperr_loclin_sameX_unif_by(SEXP XSEXP, SEXP ySEXP, SEXP gSEXP, SEXP HSEXP, SEXP nthrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat const& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec const& >::type y(ySEXP);
    Rcpp::traits::input_parameter< IntegerVector const& >::type g(gSEXP);
    Rcpp::traits::input_parameter< double const& >::type H(HSEXP);
    Rcpp::traits::input_parameter< int const >::type nthr(nthrSEXP);
    rcpp_result_gen = Rcpp::wrap(loclin_sameX_unif_by(X, y, g, H, nthr));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_looperr_fastols", (DL_FUNC) &_looperr_fastols, 4},
    {"_looperr_fastolswt", (DL_FUNC) &_looperr_fastolswt, 5},
    {"_looperr_fastols_by", (DL_FUNC) &_looperr_fastols_by, 6},
    {"_looperr_fastols_bywt", (DL_FUNC) &_looperr_fastols_bywt, 7},
    {"_looperr_loclin_diffX", (DL_FUNC) &_looperr_loclin_diffX, 6},
    {"_looperr_loclin_sameX", (DL_FUNC) &_looperr_loclin_sameX, 5},
    {"_looperr_loclin_sameX_unif", (DL_FUNC) &_looperr_loclin_sameX_unif, 3},
    {"_looperr_loclin_diffX_unif", (DL_FUNC) &_looperr_loclin_diffX_unif, 4},
    {"_looperr_loclin_sameX_unif_by", (DL_FUNC) &_looperr_loclin_sameX_unif_by, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_looperr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
