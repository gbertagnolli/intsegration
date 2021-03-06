// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// rcpp_floyd_flow
List rcpp_floyd_flow(NumericMatrix C);
RcppExport SEXP _intsegration_rcpp_floyd_flow(SEXP CSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type C(CSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_floyd_flow(C));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_floyd_flow_length
List rcpp_floyd_flow_length(NumericMatrix C, NumericMatrix L);
RcppExport SEXP _intsegration_rcpp_floyd_flow_length(SEXP CSEXP, SEXP LSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type C(CSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type L(LSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_floyd_flow_length(C, L));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_intsegration_rcpp_floyd_flow", (DL_FUNC) &_intsegration_rcpp_floyd_flow, 1},
    {"_intsegration_rcpp_floyd_flow_length", (DL_FUNC) &_intsegration_rcpp_floyd_flow_length, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_intsegration(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
