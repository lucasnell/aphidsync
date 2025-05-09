// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// make_known_fit_ptr
SEXP make_known_fit_ptr(const double& K, const arma::mat& L, const arma::vec& obs, const arma::uvec& time, const double& max_shape, const double& N0, const bool& compare_N, const double& trans_base);
RcppExport SEXP _aphidsync_make_known_fit_ptr(SEXP KSEXP, SEXP LSEXP, SEXP obsSEXP, SEXP timeSEXP, SEXP max_shapeSEXP, SEXP N0SEXP, SEXP compare_NSEXP, SEXP trans_baseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type K(KSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type obs(obsSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type time(timeSEXP);
    Rcpp::traits::input_parameter< const double& >::type max_shape(max_shapeSEXP);
    Rcpp::traits::input_parameter< const double& >::type N0(N0SEXP);
    Rcpp::traits::input_parameter< const bool& >::type compare_N(compare_NSEXP);
    Rcpp::traits::input_parameter< const double& >::type trans_base(trans_baseSEXP);
    rcpp_result_gen = Rcpp::wrap(make_known_fit_ptr(K, L, obs, time, max_shape, N0, compare_N, trans_base));
    return rcpp_result_gen;
END_RCPP
}
// known_fit_aphids0
double known_fit_aphids0(const arma::vec& pars, SEXP fit_info_ptr);
RcppExport SEXP _aphidsync_known_fit_aphids0(SEXP parsSEXP, SEXP fit_info_ptrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type fit_info_ptr(fit_info_ptrSEXP);
    rcpp_result_gen = Rcpp::wrap(known_fit_aphids0(pars, fit_info_ptr));
    return rcpp_result_gen;
END_RCPP
}
// calc_lambda
double calc_lambda(const arma::mat& L);
RcppExport SEXP _aphidsync_calc_lambda(SEXP LSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type L(LSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_lambda(L));
    return rcpp_result_gen;
END_RCPP
}
// width99
NumericVector width99(NumericVector shape);
RcppExport SEXP _aphidsync_width99(SEXP shapeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type shape(shapeSEXP);
    rcpp_result_gen = Rcpp::wrap(width99(shape));
    return rcpp_result_gen;
END_RCPP
}
// med_age
NumericVector med_age(NumericVector shape, NumericVector offset, int n_stages);
RcppExport SEXP _aphidsync_med_age(SEXP shapeSEXP, SEXP offsetSEXP, SEXP n_stagesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type shape(shapeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< int >::type n_stages(n_stagesSEXP);
    rcpp_result_gen = Rcpp::wrap(med_age(shape, offset, n_stages));
    return rcpp_result_gen;
END_RCPP
}
// logit
NumericVector logit(NumericVector p);
RcppExport SEXP _aphidsync_logit(SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(logit(p));
    return rcpp_result_gen;
END_RCPP
}
// inv_logit
NumericVector inv_logit(NumericVector a);
RcppExport SEXP _aphidsync_inv_logit(SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(inv_logit(a));
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
arma::mat make_L1(const double& shape, const double& scale, const uint32_t& n_stages, const uint32_t& adult_stage);
RcppExport SEXP _aphidsync_make_L1(SEXP shapeSEXP, SEXP scaleSEXP, SEXP n_stagesSEXP, SEXP adult_stageSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type shape(shapeSEXP);
    Rcpp::traits::input_parameter< const double& >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< const uint32_t& >::type n_stages(n_stagesSEXP);
    Rcpp::traits::input_parameter< const uint32_t& >::type adult_stage(adult_stageSEXP);
    rcpp_result_gen = Rcpp::wrap(make_L1(shape, scale, n_stages, adult_stage));
    return rcpp_result_gen;
END_RCPP
}
// sim_re
arma::vec sim_re(const arma::vec& aphids0, const arma::mat& L, const arma::uvec& time, const double& K);
RcppExport SEXP _aphidsync_sim_re(SEXP aphids0SEXP, SEXP LSEXP, SEXP timeSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type aphids0(aphids0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type time(timeSEXP);
    Rcpp::traits::input_parameter< const double& >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_re(aphids0, L, time, K));
    return rcpp_result_gen;
END_RCPP
}
// sim_N
arma::vec sim_N(const arma::vec& aphids0, const arma::mat& L, const arma::uvec& time, const double& K);
RcppExport SEXP _aphidsync_sim_N(SEXP aphids0SEXP, SEXP LSEXP, SEXP timeSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type aphids0(aphids0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type time(timeSEXP);
    Rcpp::traits::input_parameter< const double& >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_N(aphids0, L, time, K));
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
    {"_aphidsync_make_known_fit_ptr", (DL_FUNC) &_aphidsync_make_known_fit_ptr, 8},
    {"_aphidsync_known_fit_aphids0", (DL_FUNC) &_aphidsync_known_fit_aphids0, 2},
    {"_aphidsync_calc_lambda", (DL_FUNC) &_aphidsync_calc_lambda, 1},
    {"_aphidsync_width99", (DL_FUNC) &_aphidsync_width99, 1},
    {"_aphidsync_med_age", (DL_FUNC) &_aphidsync_med_age, 3},
    {"_aphidsync_logit", (DL_FUNC) &_aphidsync_logit, 1},
    {"_aphidsync_inv_logit", (DL_FUNC) &_aphidsync_inv_logit, 1},
    {"_aphidsync_beta_starts", (DL_FUNC) &_aphidsync_beta_starts, 4},
    {"_aphidsync_make_L1", (DL_FUNC) &_aphidsync_make_L1, 4},
    {"_aphidsync_sim_re", (DL_FUNC) &_aphidsync_sim_re, 4},
    {"_aphidsync_sim_N", (DL_FUNC) &_aphidsync_sim_N, 4},
    {"_aphidsync_unknown_fit_aphids0", (DL_FUNC) &_aphidsync_unknown_fit_aphids0, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_aphidsync(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
