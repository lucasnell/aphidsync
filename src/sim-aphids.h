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




#endif
