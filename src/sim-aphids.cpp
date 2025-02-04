
#include <RcppArmadillo.h>
#include <nloptrAPI.h>
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





double f(unsigned n, const double *x, double *grad, void *f_data) {

    Rcpp::checkUserInterrupt();

    OptimInfo* info = (OptimInfo*) f_data;

    arma::mat& L(info->L);
    L(0, arma::span(ADULT_STAGE-1, N_STAGES-1)) = info->fecunds * x[0];

    arma::cx_vec& eigval(info->eigval);
    arma::eig_gen(eigval, L);

    double lambda2 = arma::max(arma::abs(eigval));

    double abs_diff = std::abs(lambda2 - info->lambda);

    return abs_diff;
}







// Same as above, but for use in R:
//' @export
//' @noRd
//[[Rcpp::export]]
arma::mat make_L1(const double& shape, const double& scale) {
    arma::mat L;
    make_L1_cpp(L, shape, scale);
    return L;
}



std::vector<double> optim_L_cpp(const arma::mat& L,
                              const double& lambda,
                              const std::string& method,
                              const double& upper_bound,
                              const double& tol,
                              const int& max_iters) {

    OptimInfo info_obj(L, lambda);

    nlopt_opt opt;
    if (method == "BOBYQA") {
        opt = nlopt_create(NLOPT_LN_BOBYQA, 1);
    } else if (method == "COBYLA") {
        opt = nlopt_create(NLOPT_LN_COBYLA, 1);
    } else if (method == "NELDERMEAD") {
        opt = nlopt_create(NLOPT_LN_NELDERMEAD, 1);
    } else if (method == "SBPLX") {
        opt = nlopt_create(NLOPT_LN_SBPLX, 1);
        /*
         DIRECT_L below is a global optimizer that works particularly
         well for accuracy but in my testing takes ~14x as long as BOBYQA
         or COBYLA.
         These are probably excessive, but DIRECT_L seems to work particularly
         well. (It also takes ~14x as long as BOBYQA or COBYLA)
         */
    } else if (method == "DIRECT_L") {
        opt = nlopt_create(NLOPT_GN_DIRECT_L, 1);
    } else stop("Unsupported or wrong method\n");

    double lower_bound = 0;

    nlopt_set_lower_bounds1(opt, lower_bound);
    nlopt_set_upper_bounds1(opt, upper_bound);
    nlopt_set_min_objective(opt, f, &info_obj);
    nlopt_set_stopval(opt, 0);
    nlopt_set_ftol_rel(opt, 0);
    nlopt_set_ftol_abs(opt, 0);
    nlopt_set_xtol_rel(opt, tol);
    nlopt_set_maxeval(opt, max_iters);

    std::vector<double> x = { upper_bound / 2 };
    double opt_f;

    nlopt_result status = nlopt_optimize(opt, &(x[0]), &opt_f);

    nlopt_destroy(opt);

    x.reserve(3);
    x.push_back(opt_f);
    x.push_back(static_cast<double>(status));

    return x;

}



//' @export
//' @noRd
//[[Rcpp::export]]
std::vector<double> optim_L(const double& shape,
                            const double& scale,
                            const double& lambda,
                            const std::string& method = "BOBYQA",
                            const double& upper_bound = 1000,
                            const double& tol = 1e-4,
                            const int& max_iters = 1000) {
    arma::mat L;
    make_L1_cpp(L, shape, scale);
    std::vector<double> outs = optim_L_cpp(L, lambda, method, upper_bound,
                                           tol, max_iters);
    return outs;
}







//' @export
//' @noRd
//[[Rcpp::export]]
SEXP make_regr_ptr() {

    Environment env("package:aphidsync");
    List pf = env["poly_fits"];
    const double b0 = pf["b0"];
    const arma::vec b_shape = as<arma::vec>(pf["b_shape"]);
    const arma::vec b_scale = as<arma::vec>(pf["b_scale"]);
    const arma::vec b_lambda = as<arma::vec>(pf["b_lambda"]);
    const std::vector<arma::mat> two_way = as<std::vector<arma::mat>>(pf["two_way"]);
    const arma::cube three_way = as<arma::cube>(pf["three_way"]);

    XPtr<RegrInfo> regr_xptr(new RegrInfo(b0, b_shape, b_scale, b_lambda,
                                         two_way, three_way));

    return regr_xptr;
}


// Calculate one value of `x` (value to achieve a given lambda) from the regression
// `T` can be `arma::vec` or `subview_col<mat::elem_type>`
template <class T>
double calc_L_one__(const T& pow_shape,
                    const T& pow_scale,
                    const T& pow_lambda,
                    const double& b0,
                    const arma::vec& b_shape,
                    const arma::vec& b_scale,
                    const arma::vec& b_lambda,
                    const std::vector<arma::mat>& two_way,
                    const arma::cube& three_way) {

    double y = b0;

    for (uint32_t i = 0; i < b_shape.n_elem; i++) {

        const double& shape_i(pow_shape(i));
        y += shape_i * b_shape(i);

        for (uint32_t j = 0; j < b_scale.n_elem; j++) {

            const double& scale_j(pow_scale(j));
            if (i == 0) y += scale_j * b_scale(j);
            y += shape_i * scale_j * two_way[0](i,j);

            for (uint32_t k = 0; k < b_lambda.n_elem; k++) {

                const double& lambda_k(pow_lambda(k));
                if (j == 0) y += shape_i * lambda_k * two_way[1](i,k);
                if (i == 0) {
                    y += scale_j * lambda_k * two_way[2](j,k);
                    if (j == 0) y += lambda_k * b_lambda(k);
                }
                y += shape_i * scale_j * lambda_k * three_way(i,j,k);

            }
        }
    }

    y = std::exp(y);

    return y;

}




// Calculate Leslie matrix to reach a given lambda (using a previous polynomial
// regression).
double calc_L_one_cpp(const double& shape,
                      const double& scale,
                      const double& lambda,
                      XPtr<RegrInfo> regr_xptr) {

    const double& b0(regr_xptr->b0);
    const arma::vec& b_shape(regr_xptr->b_shape);
    const arma::vec& b_scale(regr_xptr->b_scale);
    const arma::vec& b_lambda(regr_xptr->b_lambda);
    const std::vector<arma::mat>& two_way(regr_xptr->two_way);
    const arma::cube& three_way(regr_xptr->three_way);

    arma::vec pow_shape = make_pow_vec(shape, b_shape.n_elem);
    arma::vec pow_scale = make_pow_vec(scale, b_scale.n_elem);
    arma::vec pow_lambda = make_pow_vec(lambda, b_lambda.n_elem);

    double x = calc_L_one__<arma::vec>(pow_shape, pow_scale, pow_lambda,
                                       b0, b_shape, b_scale, b_lambda,
                                       two_way, three_way);

    return x;

}






// Version outputting a vector for use in R for testing:
//' @export
//' @noRd
//[[Rcpp::export]]
arma::vec calc_L(const arma::vec& shape,
                 const arma::vec& scale,
                 const arma::vec& lambda) {

    Environment env("package:aphidsync");
    List pf = env["poly_fits"];

    const double b0 = pf["b0"];
    const arma::vec b_shape = as<arma::vec>(pf["b_shape"]);
    const arma::vec b_scale = as<arma::vec>(pf["b_scale"]);
    const arma::vec b_lambda = as<arma::vec>(pf["b_lambda"]);
    const std::vector<arma::mat> two_way = as<std::vector<arma::mat>>(pf["two_way"]);
    const arma::cube three_way = as<arma::cube>(pf["three_way"]);

    // Inputs cannot exceed bounds of the input data used for regression
    if (arma::any(shape > 5))
        stop("for this method, no shape value can be > 5");
    if (arma::any(scale > 20))
        stop("for this method, no scale value can be > 20");
    if (arma::any(lambda > 1.6))
        stop("for this method, no lambda value can be > 1.6");
    if (arma::any(lambda < 1.05))
        stop("for this method, no lambda value can be < 1.05");

    uint32_t nt = shape.n_elem;
    if (scale.n_elem != nt) stop("scale.n_elem != shape.n_elem");
    if (lambda.n_elem != nt) stop("lambda.n_elem != shape.n_elem");

    uint32_t k_shape = b_shape.n_elem;
    uint32_t k_scale = b_scale.n_elem;
    uint32_t k_lambda = b_lambda.n_elem;

    if (two_way.size() != 3)
        stop("two_way.size() != 3");
    if (two_way[0].n_rows != k_shape || two_way[0].n_cols != k_scale)
        stop("two_way[0].n_rows != k_shape || two_way[0].n_cols != k_scale");
    if (two_way[1].n_rows != k_shape || two_way[1].n_cols != k_lambda)
        stop("two_way[1].n_rows != k_shape || two_way[1].n_cols != k_lambda");
    if (two_way[2].n_rows != k_scale || two_way[2].n_cols != k_lambda)
        stop("two_way[2].n_rows != k_scale || two_way[2].n_cols != k_lambda");
    if (three_way.n_rows != k_shape || three_way.n_cols != k_scale)
        stop("three_way.n_rows != k_shape || three_way.n_cols != k_scale");
    if (three_way.n_slices != k_lambda) stop("three_way.n_slices != k_lambda");

    // Calculate powers now (they'll be used multiple times):
    arma::mat pow_shape = make_pow_mat(shape, k_shape);
    arma::mat pow_scale = make_pow_mat(scale, k_scale);
    arma::mat pow_lambda = make_pow_mat(lambda, k_lambda);

    arma::vec x(nt, arma::fill::none);

    for (uint32_t t = 0; t < nt; t++) {

        x(t) = calc_L_one__<arma::subview_col<arma::mat::elem_type>>(
            pow_shape.col(t), pow_scale.col(t), pow_lambda.col(t),
            b0, b_shape, b_scale, b_lambda, two_way, three_way);

    }

    return x;

}







//' @export
//' @noRd
//[[Rcpp::export]]
double fit_aphids0(const arma::vec& pars,
                   const arma::mat& known_L_mat,
                   const arma::vec& re,
                   const arma::uvec& time,
                   const double& fecund,
                   SEXP regr_ptr,
                   const bool& match_lambda = false,
                   const double& max_shape = 800,
                   const std::string& L_method = "BOBYQA",
                   const double& L_upper_bound = 1000,
                   const double& L_tol = 1e-4,
                   const int& L_max_iters = 1000) {

    if (arma::any(pars < 0) || pars(1) > 1) return BIG_RETURN;

    const double& shape(pars(0));
    const double& offset(pars(1));
    if (shape > max_shape) return BIG_RETURN;

    bool known_L = known_L_mat.n_rows == N_STAGES || known_L_mat.n_cols == N_STAGES;

    arma::vec aphids0 = beta_starts_cpp(shape, offset, STARTING_ABUNDANCE, N_STAGES);

    arma::mat L;
    XPtr<RegrInfo> regr_xptr(regr_ptr);
    std::vector<double> op;
    op.reserve(3);

    if (known_L) {
        L = known_L_mat;
    } else {
        if (pars.n_elem != 5) stop("pars.n_elem != 5 without known_L_mat");
        if (arma::any(pars.subvec(2, 3) == 0)) return BIG_RETURN;
        if (pars(4) > 1) return BIG_RETURN;
        // parameters for Weibull distribution for Leslie matrix:
        const double& w_shape(pars(2));
        const double& w_scale(pars(3));
        // Setup starting leslie matrix with all survivals = 1 and with
        // fecundities summing to one (lambda is also 1):
        make_L1_cpp(L, w_shape, w_scale);

        // Now use parameter #5 to estimate survivals:
        double ns;
        for (uint32_t i = 0; i < (N_STAGES-1U); i++) {
            ns = -1.0 * static_cast<double>(N_STAGES - 2U - i);
            L.diag(-1)(i) = 1 - pars(4) / std::sqrt(1 + ns * ns);
        }

        // Then either fit maximum fecundity or a given lambda (`fecund` is
        // either the max [if `!match_lambda`] or lambda [if `match_lambda`]):
        int status = 0;
        if (!match_lambda) {

            const double& max_fecund(fecund);
            double x = max_fecund;
            x /= arma::max(L.row(0).tail(N_STAGES - ADULT_STAGE + 1));
            L.row(0).tail(N_STAGES - ADULT_STAGE + 1) *= x;

        } else {

             // Adjusting lambda is less straightforward:

            const double& lambda(fecund);
            // If inputs exceed bounds of regression used to fit L, then fit
            // using an optimizer (this way is ~700x slower)
            if (w_shape > 5 || w_scale > 20 || lambda > 1.6 || lambda < 1.05) {

                op = optim_L_cpp(L, lambda, L_method, L_upper_bound,
                                 L_tol, L_max_iters);
                // Fill status and only adjust `L` if optim was successful:
                status = op[2];
                if (status >= 0) L.row(0).tail(N_STAGES - ADULT_STAGE + 1) *= op[0];

            } else {

                double x = calc_L_one_cpp(w_shape, w_scale, lambda, regr_xptr);
                L.row(0).tail(N_STAGES - ADULT_STAGE + 1) *= x;

            }
        }

        if (status < 0) return arma::datum::nan;
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



