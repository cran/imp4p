// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// fast_apply_nb_na
NumericVector fast_apply_nb_na(NumericMatrix X, int dim);
RcppExport SEXP imp4p_fast_apply_nb_na(SEXP XSEXP, SEXP dimSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type dim(dimSEXP);
    rcpp_result_gen = Rcpp::wrap(fast_apply_nb_na(X, dim));
    return rcpp_result_gen;
END_RCPP
}
// fast_apply_nb_not_na
NumericVector fast_apply_nb_not_na(NumericMatrix X, int dim);
RcppExport SEXP imp4p_fast_apply_nb_not_na(SEXP XSEXP, SEXP dimSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type dim(dimSEXP);
    rcpp_result_gen = Rcpp::wrap(fast_apply_nb_not_na(X, dim));
    return rcpp_result_gen;
END_RCPP
}
// fast_apply_sd_na_rm_T
NumericVector fast_apply_sd_na_rm_T(NumericMatrix X, int dim);
RcppExport SEXP imp4p_fast_apply_sd_na_rm_T(SEXP XSEXP, SEXP dimSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type dim(dimSEXP);
    rcpp_result_gen = Rcpp::wrap(fast_apply_sd_na_rm_T(X, dim));
    return rcpp_result_gen;
END_RCPP
}
// fast_apply_sum_na_rm_T
NumericVector fast_apply_sum_na_rm_T(NumericMatrix X, int dim);
RcppExport SEXP imp4p_fast_apply_sum_na_rm_T(SEXP XSEXP, SEXP dimSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type dim(dimSEXP);
    rcpp_result_gen = Rcpp::wrap(fast_apply_sum_na_rm_T(X, dim));
    return rcpp_result_gen;
END_RCPP
}
// fast_sim
NumericVector fast_sim(NumericVector prot, NumericMatrix mat);
RcppExport SEXP imp4p_fast_sim(SEXP protSEXP, SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type prot(protSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(fast_sim(prot, mat));
    return rcpp_result_gen;
END_RCPP
}
