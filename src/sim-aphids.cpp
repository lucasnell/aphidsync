
#include <RcppArmadillo.h>
#include <vector>
#include <math.h>
#include <algorithm>

#include "sim-aphids.h"

using namespace Rcpp;





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
arma::mat make_L1(const double& shape,
                  const double& scale,
                  const uint32_t& n_stages = 29,
                  const uint32_t& adult_stage = 9) {
    arma::mat L;
    make_L1_cpp(L, shape, scale, n_stages, adult_stage);

    return L;
}




//' Simulate per-capita growth.
//'
//' @export
//' @noRd
//[[Rcpp::export]]
arma::vec sim_re(const arma::vec& aphids0,
                 const arma::mat& L,
                 const arma::uvec& time,
                 const double& K) {
    if (aphids0.n_elem != L.n_rows || aphids0.n_elem!= L.n_cols) {
        stop("aphids0 length must equal nrows and ncols of L");
    }
    if (K <= 0) stop("K must be > 0");
    arma::vec out = sim_re_cpp(aphids0, L, time, K);
    return out;
}

//' Simulate abundance.
//'
//' @export
//' @noRd
//[[Rcpp::export]]
arma::vec sim_N(const arma::vec& aphids0,
                const arma::mat& L,
                const arma::uvec& time,
                const double& K) {
    if (aphids0.n_elem != L.n_rows || aphids0.n_elem!= L.n_cols) {
        stop("aphids0 length must equal nrows and ncols of L");
    }
    if (K <= 0) stop("K must be > 0");
    arma::vec out = sim_N_cpp(aphids0, L, time, K);
    return out;
}
