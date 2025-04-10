
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
    arma::vec pars_conv; // inverse-transformed parameters


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
        pars_conv(),
        L0(L_),
        survs(L_.diag(-1).n_elem),
        n_pars(0) {};

    void adjust_L() {
        L = L0;
        L *= pars_conv(2);
        return;
    }
    void adjust_fecunds() {
        L = L0;
        L.row(0) *= pars_conv(2);
        return;
    }
    void adjust_survs() {
        L = L0;
        survs = L.diag(-1);
        logit__(survs);
        survs += pars_conv(3);
        inv_logit__(survs);
        L.diag(-1) = survs;
        return;
    }

    /*
     Convert `pars` from scale of fitting (used in optimizer) to scale
     of interest (used `sim_re_cpp` and `sim_N_cpp`):
     Note: This does not transform the 4th parameter (if there is one) bc
     it can take any value.
     See description above `known_fit_aphids0` or the inside of function
     `KnownLeslieFitInfo::adjust_survs` for more info.
     */
    void convert_pars(const arma::vec& pars) {

        if (n_pars != pars.n_elem) {
            n_pars = pars.n_elem;
            if (n_pars != 2 && n_pars != 3 && n_pars != 4) {
                stop("length(pars) should be 2, 3, or 4");
            }
            pars_conv.set_size(pars.n_elem);
        }

        // This keeps shapes above 1:
        pars_conv(0) = std::pow(trans_base, pars(0)) + 1.0;
        // This keeps offset between 0 and 1:
        pars_conv(1) = inv_logit__(pars(1));
        // This keeps Leslie or fecundity multiple above 0:
        if (n_pars > 2) pars_conv(2) = std::pow(trans_base, pars(2));
        if (n_pars > 3) pars_conv(3) = pars(3);

        return;
    }



private:

    arma::mat L0;
    arma::vec survs;
    uint32_t n_pars;

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
//' 1. shape = trans_base^(pars(0)) + 1
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

    // Convert from scale used by optimizer to scale used here
    // (stored in `fit_info.pars_conv`):
    fit_info.convert_pars(pars);
    const double& shape(fit_info.pars_conv(0));
    const double& max_shape(fit_info.max_shape);
    if (shape > max_shape) return BIG_RETURN;
    const double& offset(fit_info.pars_conv(1));

    // Scale entire Leslie matrix if 3 parameters provide:
    if (pars.n_elem == 3) {
        fit_info.adjust_L(); // will use value in `fit_info.pars_conv`
    }
    // Scale fecundities and survivals separately if 4 parameters provide:
    if (pars.n_elem == 4) {
        // both of these will use values in `fit_info.pars_conv`
        fit_info.adjust_fecunds();
        fit_info.adjust_survs();
    }

    const double& K(fit_info.K);
    const arma::mat& L(fit_info.L);
    const arma::vec& obs(fit_info.obs);
    const arma::uvec& time(fit_info.time);
    const double& N0(fit_info.N0);
    const bool& compare_N(fit_info.compare_N);

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
