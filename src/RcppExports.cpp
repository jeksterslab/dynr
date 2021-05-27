// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;


// rcpp_saem_interface
Rcpp::List rcpp_saem_interface(Rcpp::List model_sexp, Rcpp::List data_sexp, bool weight_flag, bool debug_flag, bool optimization_flag, bool hessian_flag, bool verbose);


RcppExport SEXP wrap_rcpp_saem_interface(SEXP model_sexpSEXP, SEXP data_sexpSEXP, SEXP weight_flagSEXP, SEXP debug_flagSEXP, SEXP optimization_flagSEXP, SEXP hessian_flagSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type model_sexp(model_sexpSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type data_sexp(data_sexpSEXP);
    Rcpp::traits::input_parameter< bool >::type weight_flag(weight_flagSEXP);
    Rcpp::traits::input_parameter< bool >::type debug_flag(debug_flagSEXP);
    Rcpp::traits::input_parameter< bool >::type optimization_flag(optimization_flagSEXP);
    Rcpp::traits::input_parameter< bool >::type hessian_flag(hessian_flagSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_saem_interface(model_sexp, data_sexp, weight_flag, debug_flag, optimization_flag, hessian_flag, verbose));
    return rcpp_result_gen;
END_RCPP
}


