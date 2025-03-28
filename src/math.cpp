

#include <RcppArmadillo.h>
#include <vector>
#include <math.h>
#include <algorithm>
#include "math.h"
#include "sim-aphids.h"



using namespace Rcpp;



//' Calculate lambda (exponential growth rate) from Leslie matrix.
//'
//' @param L Leslie matrix.
//'
//' @export
//'
//[[Rcpp::export]]
double calc_lambda(const arma::mat& L) {
    if (L.n_rows != L.n_cols) stop("Leslie matrix must be square");
    double lambda = calc_lambda_cpp(L);
    return lambda;
}


/*
 =====================================================================================
 =====================================================================================
  Summary stats
 =====================================================================================
 =====================================================================================
 */

//' Width of 99th quartile
//'
//' @param shape Vector of shape parameters for the symmetrical beta
//'     distribution of abundances.
//'
//' @export
//'
//[[Rcpp::export]]
NumericVector width99(NumericVector shape) {
    NumericVector out(shape.size());
    double lo, hi;
    for (uint32_t i = 0; i < shape.size(); i++) {
        lo = R::qbeta(0.005, shape(i), shape(i), true, false);
        hi = R::qbeta(0.995, shape(i), shape(i), true, false);
        out(i) = std::abs(hi - lo);
    }
    return out;
}


//' Median age
//'
//' @inheritParams width99
//' @param offset Vector of offset value(s) for the symmetrical beta
//'     distribution of abundances.
//'     Must be the same length as `shape`.
//'
//' @export
//'
//[[Rcpp::export]]
NumericVector med_age(NumericVector shape,
                      NumericVector offset,
                      int n_stages = 29) {
    uint32_t n = shape.size();
    if (offset.size() != n) stop("shape and offset must be same length");
    if (n_stages < 2) stop("n_stages must be > 1");
    NumericVector out(n);
    arma::vec abunds(n_stages, arma::fill::none);
    double idx;
    for (uint32_t i = 0; i < n; i++) {
        abunds = beta_starts_cpp(shape(i), offset(i), 1, n_stages);
        uint32_t j = 0;
        while (abunds(j) <= 0.5) {
            j++;
            abunds(j) += abunds(j-1);
        }
        idx = static_cast<double>(j);
        if (j > 0 && abunds(j-1) == 0.5) {
            // This is for rare occurrence where median age is exactly between
            // two age groups:
            idx += static_cast<double>(j-1);
            idx /= 2.0;
        }
        out(i) = idx + 1.0;
    }
    return out;
}



/*
 =====================================================================================
 =====================================================================================
 Logit and inverse logit
 =====================================================================================
 =====================================================================================
 */

//' Logit and inverse logit functions.
//'
//'
//' @name logit
//' @export
//'
//[[Rcpp::export]]
NumericVector logit(NumericVector p) {
    NumericVector out(p.size());
    for (uint32_t i = 0; i < p.size(); i++) {
        logit__(p[i], out[i]);
    }
    return out;
}
//' @describeIn logit
//'
//' @export
//'
//[[Rcpp::export]]
NumericVector inv_logit(NumericVector a){
    NumericVector out(a.size());
    for (uint32_t i = 0; i < a.size(); i++) {
        inv_logit__(a[i], out[i]);
    }
    return out;
}
