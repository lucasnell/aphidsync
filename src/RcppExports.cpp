// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// known_fit_aphids0
double known_fit_aphids0(const arma::vec& pars, const arma::mat& L, const arma::vec& obs, const arma::uvec& time, const double& max_shape, const bool& compare_N);
RcppExport SEXP _aphidsync_known_fit_aphids0(SEXP parsSEXP, SEXP LSEXP, SEXP obsSEXP, SEXP timeSEXP, SEXP max_shapeSEXP, SEXP compare_NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type time(timeSEXP);
    Rcpp::traits::input_parameter< const double& >::type max_shape(max_shapeSEXP);
    Rcpp::traits::input_parameter< const bool& >::type compare_N(compare_NSEXP);
    rcpp_result_gen = Rcpp::wrap(known_fit_aphids0(pars, L, obs, time, max_shape, compare_N));
    return rcpp_result_gen;
END_RCPP
}
// beta_starts
arma::vec beta_starts(const double& shape, const double& offset, const uint32_t& total_aphids0, const uint32_t& compartments);
RcppExport SEXP _aphidsync_beta_starts(SEXP shapeSEXP, SEXP offsetSEXP, SEXP total_aphids0SEXP, SEXP compartmentsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type shape(shapeSEXP);
    Rcpp::traits::input_parameter< const double& >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< const uint32_t& >::type total_aphids0(total_aphids0SEXP);
    Rcpp::traits::input_parameter< const uint32_t& >::type compartments(compartmentsSEXP);
    rcpp_result_gen = Rcpp::wrap(beta_starts(shape, offset, total_aphids0, compartments));
    return rcpp_result_gen;
END_RCPP
}
// make_L1
arma::mat make_L1(const double& shape, const double& scale);
RcppExport SEXP _aphidsync_make_L1(SEXP shapeSEXP, SEXP scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type shape(shapeSEXP);
    Rcpp::traits::input_parameter< const double& >::type scale(scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(make_L1(shape, scale));
    return rcpp_result_gen;
END_RCPP
}
// unknown_fit_aphids0
double unknown_fit_aphids0(const arma::vec& pars, const arma::vec& re, const arma::uvec& time, const double& max_f, const bool& fit_survs, const double& max_shape);
RcppExport SEXP _aphidsync_unknown_fit_aphids0(SEXP parsSEXP, SEXP reSEXP, SEXP timeSEXP, SEXP max_fSEXP, SEXP fit_survsSEXP, SEXP max_shapeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type re(reSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type time(timeSEXP);
    Rcpp::traits::input_parameter< const double& >::type max_f(max_fSEXP);
    Rcpp::traits::input_parameter< const bool& >::type fit_survs(fit_survsSEXP);
    Rcpp::traits::input_parameter< const double& >::type max_shape(max_shapeSEXP);
    rcpp_result_gen = Rcpp::wrap(unknown_fit_aphids0(pars, re, time, max_f, fit_survs, max_shape));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_aphidsync_known_fit_aphids0", (DL_FUNC) &_aphidsync_known_fit_aphids0, 6},
    {"_aphidsync_beta_starts", (DL_FUNC) &_aphidsync_beta_starts, 4},
    {"_aphidsync_make_L1", (DL_FUNC) &_aphidsync_make_L1, 2},
    {"_aphidsync_unknown_fit_aphids0", (DL_FUNC) &_aphidsync_unknown_fit_aphids0, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_aphidsync(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
