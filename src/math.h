# ifndef __APHIDSYNC_MATH_H
# define __APHIDSYNC_MATH_H


#include <RcppArmadillo.h>
#include <vector>
#include <math.h>
#include <algorithm>


using namespace Rcpp;





// Calculate lambda from Leslie matrix:
inline double calc_lambda_cpp(const arma::mat& L) {

    arma::cx_vec eigval;
    arma::eig_gen(eigval, L);

    double lambda = arma::max(arma::abs(eigval));

    return lambda;

}


/*
 =====================================================================================
 =====================================================================================
 Logit and inverse logit
 =====================================================================================
 =====================================================================================
 */

inline void logit__(const double& p, double& out) {
    out = std::log(p / (1-p));
    return;
}
inline void inv_logit__(const double& a, double& out) {
    out = 1 / (1 + std::exp(-a));
    return;
}
inline double logit__(const double& p) {
    double out = std::log(p / (1-p));
    return out;
}
inline double inv_logit__(const double& a) {
    double out = 1 / (1 + std::exp(-a));
    return out;
}
// Overloaded for editing arma vectors in place
inline void logit__(arma::vec& p_vec) {
    double out;
    for (double& p : p_vec) {
        out = std::log(p / (1-p));
        p = out;
    }
    return;
}
inline void inv_logit__(arma::vec& a_vec) {
    double out;
    for (double& a : a_vec) {
        out = 1 / (1 + std::exp(-a));
        a = out;
    }
    return;
}


#endif



