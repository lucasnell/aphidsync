
#include <RcppArmadillo.h>
#include <vector>
#include <math.h>
#include <algorithm>

#include "sim-aphids.h"
#include "math.h"

using namespace Rcpp;


struct KnownLeslieFitInfo {

    double K;
    arma::mat L;
    arma::vec obs;
    arma::uvec time;
    double max_shape;
    double N0;
    bool compare_N;
    double trans_base;


    KnownLeslieFitInfo(const double& K_,
                       const arma::mat& L_,
                       const arma::vec& obs_,
                       const arma::uvec& time_,
                       const double& max_shape_,
                       const double& N0_,
                       const bool& compare_N_,
                       const double& trans_base_) :
        K(K_),
        L(L_),
        obs(obs_),
        time(time_),
        max_shape(max_shape_),
        N0(N0_),
        compare_N(compare_N_),
        trans_base(trans_base_),
        L0(L_),
        survs(L_.diag(-1).n_elem) {};

    void adjust_L(const double& x) {
        L = L0;
        double z = std::pow(trans_base, x);
        L *= z;
        return;
    }
    void adjust_fecunds(const double& x) {
        L = L0;
        double z = std::pow(trans_base, x);
        L.row(0) *= z;
        return;
    }
    void adjust_survs(const double& x) {
        L = L0;
        survs = L.diag(-1);
        logit__(survs);
        survs += x;
        inv_logit__(survs);
        L.diag(-1) = survs;
        return;
    }


private:

    arma::mat L0;
    arma::vec survs;

};



//' Make object necessary for use in `known_fit_aphids0`.
//' Note: this should NOT be shared across threads, only for one use of
//' an optimization process.
//'
//' @export
//' @noRd
//[[Rcpp::export]]
SEXP make_known_fit_ptr(const double& K,
                        const arma::mat& L,
                        const arma::vec& obs,
                        const arma::uvec& time,
                        const double& max_shape,
                        const double& N0,
                        const bool& compare_N,
                        const double& trans_base) {

    if (obs.n_elem != time.n_elem) stop("obs.n_elem != time.n_elem");
    if (L.n_rows != L.n_cols) stop("L.n_rows != L.n_cols");

    XPtr<KnownLeslieFitInfo> fit_info_xptr(
            new KnownLeslieFitInfo(K, L, obs, time, max_shape, N0, compare_N,
                                   trans_base), true);

    return fit_info_xptr;


}

//' This function scales parameters as follows (in order):
//'
//' 1. shape = trans_base^(pars(0))
//' 2. offset = inv_logit(pars(1))
//' 3. fecund_x = trans_base^(pars(2))
//' 4. surv_x = inv_logit(pars(3))**
//'
//' ** surv_x is added to logit-transformed survivals before back-transforming
//'    using inv_logit, so it's effectively inv_logit transformed here.
//'
//' @export
//' @noRd
//[[Rcpp::export]]
double known_fit_aphids0(const arma::vec& pars,
                         SEXP fit_info_ptr) {

    if (pars.n_elem < 2) stop("pars.n_elem < 2 in known_fit_aphids0");
    // if (arma::any(pars < 0) || pars(1) > 1) return BIG_RETURN;

    XPtr<KnownLeslieFitInfo> fit_info_xptr(fit_info_ptr);
    KnownLeslieFitInfo& fit_info(*fit_info_xptr);

    // Scale entire Leslie matrix if 3 parameters provide:
    if (pars.n_elem == 3) {
        fit_info.adjust_L(pars(2)); // transformed inside `adjust_L`
    }
    // Scale fecundities and survivals separately if 4 parameters provide:
    if (pars.n_elem == 4) {
        fit_info.adjust_fecunds(pars(2)); // transformed inside `adjust_fecunds`
        fit_info.adjust_survs(pars(3));
    }

    const double& K(fit_info.K);
    const arma::mat& L(fit_info.L);
    const arma::vec& obs(fit_info.obs);
    const arma::uvec& time(fit_info.time);
    const double& max_shape(fit_info.max_shape);
    const double& N0(fit_info.N0);
    const bool& compare_N(fit_info.compare_N);

    double shape = std::exp(pars(0));
    if (shape > max_shape) return BIG_RETURN;
    double offset = inv_logit__(pars(1));

    uint32_t n_stages = L.n_rows;

    arma::vec aphids0 = beta_starts_cpp(shape, offset, N0, n_stages);

    // Simulate with given parameters (either per-capita growth or abundance):
    arma::vec pred;
    if (compare_N) {
        pred = sim_N_cpp(aphids0, L, time, K);
    } else pred = sim_re_cpp(aphids0, L, time, K);

    double sae = 0;
    for (uint32_t t = 0; t < pred.n_elem; t++) {
        sae += std::abs(pred(t) - obs(t));
    }
    if (sae > BIG_RETURN) sae = BIG_RETURN;

    return sae;

}
