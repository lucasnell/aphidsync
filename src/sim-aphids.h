# ifndef __APHIDSYNC_SIM_APHIDS_H
# define __APHIDSYNC_SIM_APHIDS_H


#include <RcppArmadillo.h>
#include <vector>
#include <math.h>
#include <algorithm>


using namespace Rcpp;


constexpr double BIG_RETURN = 1e10;
constexpr uint32_t N_STAGES = 29;
constexpr uint32_t ADULT_STAGE = 9;
constexpr double STARTING_ABUNDANCE = 32;
constexpr double REAL_K = 1800;








// Create Leslie matrix with only survivals = 1 and with fecundities that
// sum to 1. Lambda for this matrix is also 1.
inline void make_L1_cpp(arma::mat& L,
                        const double& shape,
                        const double& scale) {
    if (L.n_rows != N_STAGES || L.n_cols != N_STAGES)
        L.set_size(N_STAGES, N_STAGES);
    L.zeros();
    L.diag(-1).fill(1);
    double pweib_val, pweib_val0;
    double pweib_diff_sum = 0;
    uint32_t adult_stages = N_STAGES - ADULT_STAGE + 1;
    arma::rowvec pweib_diffs(adult_stages, arma::fill::none);
    for (uint32_t i = 0; i <= adult_stages; i++) {
        pweib_val = R::pweibull(i, shape, scale, true, false);
        if (i > 0) {
            pweib_diffs(i-1) = pweib_val - pweib_val0;
            pweib_diff_sum += pweib_diffs(i-1);
        }
        pweib_val0 = pweib_val;
    }
    // So they sum to one:
    pweib_diffs /= pweib_diff_sum;
    L.row(0).tail(adult_stages) = pweib_diffs;

    return;
}




// Internal C++ code for use in functions exported to R
inline arma::vec beta_starts_cpp(const double& shape,
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



// Simulate values of per-capita growth to be compared to observed values
inline arma::vec sim_re(const arma::vec& aphids0,
                        const arma::mat& L,
                        const arma::uvec& time,
                        const double& K) {

    uint32_t max_t = arma::max(time) + 1;

    arma::vec re_pred(max_t - 1, arma::fill::none);

    arma::vec Nt = aphids0;
    arma::vec dN(Nt.n_elem, arma::fill::none);
    double Nsums0 = arma::accu(aphids0);
    double S, Nsums;
    for (uint32_t t = 1; t < max_t; t++) {
        double& z(Nsums0);
        S = 1 - z / K;
        dN = S * ((L * Nt) - Nt);
        Nt += dN;
        Nsums = arma::accu(Nt);
        re_pred(t-1) = Nsums / Nsums0;
        Nsums0 = Nsums;
    }

    // If `re_pred` doesn't match length of `time` -1, then use `time` as a
    // vector of indices
    if (re_pred.n_elem != (time.n_elem-1)) {
        re_pred = re_pred(time.head(time.n_elem-1));
    }

    return re_pred;

}


// Simulate abundances (summed across stages) to be compared to observed values
inline arma::vec sim_N(const arma::vec& aphids0,
                       const arma::mat& L,
                       const arma::uvec& time,
                       const double& K) {

    uint32_t max_t = arma::max(time) + 1;

    arma::vec N_pred(max_t, arma::fill::none);

    arma::vec Nt = aphids0;
    arma::vec dN(Nt.n_elem, arma::fill::none);
    N_pred(0) = arma::accu(aphids0);
    double S;
    for (uint32_t t = 1; t < max_t; t++) {
        double& z(N_pred(t-1));
        S = 1 - z / K;
        dN = S * ((L * Nt) - Nt);
        Nt += dN;
        N_pred(t) = arma::accu(Nt);
    }

    // If `N_pred` doesn't match length of `time`, then use `time` as a
    // vector of indices
    if (N_pred.n_elem != time.n_elem) N_pred = N_pred(time);

    return N_pred;

}



#endif
