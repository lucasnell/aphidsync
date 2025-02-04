# ifndef __APHIDSYNC_SIM_APHIDS_H
# define __APHIDSYNC_SIM_APHIDS_H


#include <RcppArmadillo.h>
#include <nloptrAPI.h>
#include <vector>
#include <math.h>
#include <algorithm>


using namespace Rcpp;


constexpr double BIG_RETURN = 1e10;
constexpr uint32_t N_STAGES = 29;
constexpr uint32_t ADULT_STAGE = 9;
constexpr double STARTING_ABUNDANCE = 32;
constexpr double REAL_K = 1800;



inline double logit_cpp(const double& p) {
    return std::log(p / (1-p));
}

inline double inv_logit_cpp(const double& a) {;
    return 1 / (1 + std::exp(-a));
}




struct OptimInfo {

    arma::rowvec fecunds;
    arma::mat L;
    double lambda;
    arma::cx_vec eigval;

    OptimInfo(const arma::mat& L_, const double& lambda_)
        : fecunds(L_.row(0).subvec(ADULT_STAGE-1, N_STAGES-1)),
          L(L_),
          lambda(lambda_),
          eigval() {};

};




struct RegrInfo {
    double b0;
    arma::vec b_shape;
    arma::vec b_scale;
    arma::vec b_lambda;
    std::vector<arma::mat> two_way;
    arma::cube three_way;

    RegrInfo(const double& b0_,
             const arma::vec& b_shape_,
             const arma::vec& b_scale_,
             const arma::vec& b_lambda_,
             const std::vector<arma::mat>& two_way_,
             const arma::cube& three_way_)
        : b0(b0_),
          b_shape(b_shape_),
          b_scale(b_scale_),
          b_lambda(b_lambda_),
          two_way(two_way_),
          three_way(three_way_) {};

};




inline arma::mat make_pow_mat(const arma::vec& z, const uint32_t& mp) {

    arma::rowvec zt = z.t();
    arma::mat pm(mp, z.n_elem, arma::fill::none);
    pm.row(0) = zt;
    for (uint32_t j = 1; j < mp; j++) {
        pm.row(j) = pm.row(j-1) % zt;
    }
    return pm;

}

inline arma::vec make_pow_vec(const double& z, const uint32_t& mp) {

    arma::vec pm(mp, arma::fill::none);
    pm(0) = z;
    for (uint32_t j = 1; j < mp; j++) {
        pm(j) = pm(j-1) * z;
    }
    return pm;

}




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




#endif
