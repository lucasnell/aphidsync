
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
arma::mat make_L1(const double& shape, const double& scale) {
    arma::mat L;
    make_L1_cpp(L, shape, scale);
    return L;
}

