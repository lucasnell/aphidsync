
#include <RcppArmadillo.h>
#include <vector>
#include <math.h>
#include <algorithm>

#include "sim-aphids.h"

using namespace Rcpp;




// Internal C++ code for use in functions exported to R
arma::vec beta_starts_cpp(const double& shape,
                          const double& offset,
                          const uint32_t& total_aphids0,
                          const uint32_t& compartments) {

    double dbl_total_aphids0 = static_cast<double>(total_aphids0);

    arma::vec times = arma::linspace<arma::vec>(0, 1, compartments+1U);
    times += offset;

    arma::vec aphids0(compartments, arma::fill::none);
    double sum_aphids0 = 0;

    double corr_time, pbeta_val, pbeta_val0;
    for (uint32_t i = 0; i <= compartments; i++) {
        corr_time = times(i);
        if (corr_time > 1) corr_time -= 1;
        pbeta_val = R::pbeta(corr_time, shape, shape, true, false);
        if (times(i) > 1) pbeta_val += 1;
        if (i > 0) {
            aphids0(i-1) = dbl_total_aphids0 * (pbeta_val - pbeta_val0);
            sum_aphids0 += aphids0(i-1);
        }
        pbeta_val0 = pbeta_val;
    }

    double sum_diff = std::round(sum_aphids0) != dbl_total_aphids0;
    if (sum_diff != 0) {
        aphids0 *= (dbl_total_aphids0 / sum_aphids0);
        sum_diff = std::round(arma::accu(aphids0)) != dbl_total_aphids0;
    }
    if (sum_diff != 0) {
        std::string err = "beta_starts magnitude error (";
        err += std::to_string(sum_diff) + ", shape = ";
        err += std::to_string(shape) + ", offset = ";
        err += std::to_string(offset) + ")";
        Rcpp::warning(err.c_str());
    }
    return aphids0;
}

//' @export
//' @noRd
//[[Rcpp::export]]
arma::vec beta_starts(const double& shape,
                      const double& offset,
                      const uint32_t& total_aphids0,
                      const uint32_t& compartments = 29) {

    if (shape <= 0) stop("shape <= 0");
    if (offset < 0) stop("offset < 0");
    if (offset > 1) stop("offset > 1");
    if (compartments == 0) stop("compartments == 0");
    if (total_aphids0 == 0) stop("total_aphids0 == 0");

    arma::vec aphids0 = beta_starts_cpp(shape, offset, total_aphids0,
                                        compartments);
    return aphids0;
}








// Create Leslie matrix with only survivals = 1 and with fecundities that
// sum to 1. Lambda for this matrix is also 1. For use in R.
//' @export
//' @noRd
//[[Rcpp::export]]
arma::mat make_L1(const double& shape, const double& scale) {
    arma::mat L;
    make_L1_cpp(L, shape, scale);
    return L;
}








//' @export
//' @noRd
//[[Rcpp::export]]
double fit_aphids0(const arma::vec& pars,
                   const arma::mat& known_L_mat,
                   const arma::vec& re,
                   const arma::uvec& time,
                   const double& max_f,
                   const bool& fit_survs,
                   const double& max_shape) {

    if (arma::any(pars < 0) || pars(1) > 1) return BIG_RETURN;

    const double& shape(pars(0));
    const double& offset(pars(1));
    if (shape > max_shape) return BIG_RETURN;

    bool known_L = known_L_mat.n_rows == N_STAGES || known_L_mat.n_cols == N_STAGES;

    arma::vec aphids0 = beta_starts_cpp(shape, offset, STARTING_ABUNDANCE, N_STAGES);

    arma::mat L;

    if (known_L) {
        L = known_L_mat;
    } else {
        if (pars.n_elem < 4) stop("pars.n_elem < 4 without known_L_mat");
        if (arma::any(pars.subvec(2, 3) == 0)) return BIG_RETURN;
        if (fit_survs && pars.n_elem < 5) stop("fit_survs && pars.n_elem < 5");
        // parameters for Weibull distribution for Leslie matrix:
        const double& w_shape(pars(2));
        const double& w_scale(pars(3));
        // Setup starting leslie matrix with all survivals = 1 and with
        // fecundities summing to one (lambda is also 1):
        make_L1_cpp(L, w_shape, w_scale);

        // Now use parameter #5 to estimate survivals:
        if (fit_survs) {
            if (pars(4) > 1) return BIG_RETURN;
            double ns;
            for (uint32_t i = 0; i < (N_STAGES-1U); i++) {
                ns = -1.0 * static_cast<double>(N_STAGES - 2U - i);
                L.diag(-1)(i) = 1 - pars(4) / std::sqrt(1 + ns * ns);
            }
        }

        double x = max_f;
        x /= arma::max(L.row(0).tail(N_STAGES - ADULT_STAGE + 1));
        L.row(0).tail(N_STAGES - ADULT_STAGE + 1) *= x;

    }

    // Extrapolate through time:
    uint32_t max_t = arma::max(time) + 1;
    arma::mat Nmat(aphids0.n_elem, max_t, arma::fill::none);
    arma::vec Nsums(max_t, arma::fill::none);
    Nmat.col(0) = aphids0;
    Nsums(0) = arma::accu(aphids0);
    arma::vec re_pred(max_t - 1, arma::fill::none);
    double S;
    for (uint32_t t = 1; t < max_t; t++) {
        double& z(Nsums(t-1));
        S = 1 / (1 + z / REAL_K);
        Nmat.col(t) = S * (L * Nmat.col(t-1));
        Nsums(t) = arma::accu(Nmat.col(t));
        re_pred(t-1) = Nsums(t) / Nsums(t-1);
    }

    if (re_pred.n_elem != (time.n_elem-1)) {
        re_pred = re_pred(time.head(time.n_elem-1));
    }
    if (re.n_elem != time.n_elem) stop("re.n_elem != time.n_elem");

    double sae = 0;
    for (uint32_t t = 0; t < re_pred.n_elem; t++) {
        sae += std::abs(re_pred(t) - re(t));
    }
    if (sae > BIG_RETURN) sae = BIG_RETURN;

    return sae;

}



