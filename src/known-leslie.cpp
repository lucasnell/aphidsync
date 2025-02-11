
#include <RcppArmadillo.h>
#include <vector>
#include <math.h>
#include <algorithm>

#include "sim-aphids.h"

using namespace Rcpp;



//' @export
//' @noRd
//[[Rcpp::export]]
double known_fit_aphids0(const arma::vec& pars,
                         const arma::mat& L,
                         const arma::vec& re,
                         const arma::uvec& time,
                         const double& max_shape) {

     if (pars.n_elem < 3) stop("pars.n_elem < 3 in known_fit_aphids0");

     if (arma::any(pars < 0) || pars(1) > 1) return BIG_RETURN;

     const double& shape(pars(0));
     const double& offset(pars(1));
     const double& K(pars(2));
     if (shape > max_shape) return BIG_RETURN;

     if(L.n_rows != N_STAGES || L.n_cols != N_STAGES)
         stop("L.n_rows != N_STAGES || L.n_cols != N_STAGES");

     arma::vec aphids0 = beta_starts_cpp(shape, offset, STARTING_ABUNDANCE, N_STAGES);

     // Simulate with given parameters:
     arma::vec re_pred = sim_re(aphids0, L, time, K);

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
