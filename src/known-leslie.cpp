
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
                         const arma::vec& obs,
                         const arma::uvec& time,
                         const double& max_shape,
                         const bool& compare_N = false) {

     if (pars.n_elem < 3) stop("pars.n_elem < 3 in known_fit_aphids0");

     if (arma::any(pars < 0) || pars(1) > 1) return BIG_RETURN;

     const double& shape(pars(0));
     const double& offset(pars(1));
     const double& K(pars(2));
     if (shape > max_shape) return BIG_RETURN;

     if(L.n_rows != N_STAGES || L.n_cols != N_STAGES)
         stop("L.n_rows != N_STAGES || L.n_cols != N_STAGES");

     arma::vec aphids0 = beta_starts_cpp(shape, offset, STARTING_ABUNDANCE, N_STAGES);

     // Simulate with given parameters (either per-capita growth or abundance):
     arma::vec pred;
     if (compare_N) {
         pred = sim_N(aphids0, L, time, K);
     } else pred = sim_re(aphids0, L, time, K);

     if (pred.n_elem != (time.n_elem-1)) {
         pred = pred(time.head(time.n_elem-1));
     }
     if (obs.n_elem != time.n_elem) stop("obs.n_elem != time.n_elem");

     double sae = 0;
     for (uint32_t t = 0; t < pred.n_elem; t++) {
         sae += std::abs(pred(t) - obs(t));
     }
     if (sae > BIG_RETURN) sae = BIG_RETURN;

     return sae;

 }
