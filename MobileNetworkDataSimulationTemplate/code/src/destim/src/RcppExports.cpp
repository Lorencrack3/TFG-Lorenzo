// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// cppfunique
NumericMatrix cppfunique(NumericMatrix mat, double tol);
RcppExport SEXP _destim_cppfunique(SEXP matSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type mat(matSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(cppfunique(mat, tol));
    return rcpp_result_gen;
END_RCPP
}
// cppfuniqueind
IntegerVector cppfuniqueind(NumericMatrix mat, double tol);
RcppExport SEXP _destim_cppfuniqueind(SEXP matSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type mat(matSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(cppfuniqueind(mat, tol));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_destim_cppfunique", (DL_FUNC) &_destim_cppfunique, 2},
    {"_destim_cppfuniqueind", (DL_FUNC) &_destim_cppfuniqueind, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_destim(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
