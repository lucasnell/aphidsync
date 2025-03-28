
#include <RcppArmadillo.h>
#include <vector>
#include <math.h>
#include <algorithm>

#include "sim-aphids.h"

using namespace Rcpp;



//' @export
//' @noRd
//[[Rcpp::export]]
double unknown_fit_aphids0(const arma::vec& pars,
                           const arma::vec& re,
                           const arma::uvec& time,
                           const double& max_f,
                           const bool& fit_survs,
                           const double& max_shape) {

     if (pars.n_elem < 4) stop("pars.n_elem < 4 in unknown_fit_aphids0");

     if (arma::any(pars < 0) || pars(1) > 1) return BIG_RETURN;

     const double& shape(pars(0));
     const double& offset(pars(1));
     if (shape > max_shape) return BIG_RETURN;

     arma::vec aphids0 = beta_starts_cpp(shape, offset, STARTING_ABUNDANCE, N_STAGES);

     arma::mat L;

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

     // Simulate with given parameters:
     arma::vec re_pred = sim_re_cpp(aphids0, L, time, REAL_K);

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
