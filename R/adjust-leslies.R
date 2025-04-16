
#' Adjust Leslie matrix to match a given overall growth rate (lambda)
#'
#' @param L Square numeric Leslie matrix that will be adjusted.
#' @param lambda Single numeric indicating lambda that output matrix will
#'     be adjusted to.
#' @param adjust_part Single character indicating which part of the Leslie
#'     matrix to adjust. Options are survivals (`"s"`) or fecundities (`"f"`).
#'     Defaults to `"f"`.
#' @param max_val Single numeric indicating the maximum value of the modifier
#'     to use to achieve the desired lambda.
#'     Defaults to `NA`, which results in a very coarse search (by powers of
#'     10) for a maximum.
#' @param direct_l Single logical for whether to use the `DIRECT_L` global
#'     optimizer from `nloptr`. It's usually more accurate, but slower.
#'     Defaults to `FALSE`.
#' @param precision Single numeric indicating how close the output lambda must
#'     be to the desired value.
#'     Defaults to `1e-4`.
#' @param optim_tol Single numeric indicating the tolerance to use for the
#'     optimization process.
#'     If `direct_l = FALSE`, this is the `tol` argument for `stats::optimize`.
#'     If `direct_l = TRUE`, this is used for `xtol_rel` in NLOPT options
#'     (see `?nl.opts`).
#'     Defaults to `.Machine$double.eps^0.5`.
#'
#' @returns An adjusted Leslie matrix.
#'
#' @export
#'
#' @importFrom stats optimize
#' @importFrom nloptr directL
#'
adjust_lambda <- function(L,
                          lambda,
                          adjust_part = "f",
                          max_val = NA,
                          direct_l = FALSE,
                          precision = 1e-4,
                          optim_tol = .Machine$double.eps^0.5) {

    type_checker(L, "L", "matrix", l = NA, .min = 0)
    if (nrow(L) != ncol(L)) stop("\nERROR: L must be square")
    if (calc_lambda(L) <= 0) stop("\nERROR: L cannot start with lambda = 0")
    type_checker(lambda, "lambda", "numeric", .min = .Machine$double.eps)
    type_checker(adjust_part, "adjust_part", "character")
    type_checker(direct_l, "direct_l", "logical")
    type_checker(precision, "precision", "numeric", .min = 0)
    type_checker(optim_tol, "optim_tol", "numeric", .min = 0)

    adjust_part <- match.arg(tolower(adjust_part), c("f", "s"))

    if (adjust_part == "f") {
        new_L_fun <- function(x) {
            L[1,] <- L[1,] * x
            return(L)
        }
    } else {
        L2 <- L
        L2[row(L2) - col(L2) == 1 & L2 > 0] <- 1
        if (calc_lambda(L2) < lambda)
            stop(paste0("\nERROR: desired lambda (", lambda, ") not achievable",
                       " with only survivals. The max possible is ",
                       calc_lambda(L2)))
        new_L_fun <- function(x) {
            survs <- L[row(L) - col(L) == 1]
            survs <- inv_logit(logit(survs) + x)
            L[row(L) - col(L) == 1] <- survs
            return(L)
        }
    }

    op_fun <- function(x) {
        L <- new_L_fun(x)
        return(abs(calc_lambda(L) - lambda))
    }

    if (is.na(max_val)) {
        max_val <- 1
        while (calc_lambda(new_L_fun(max_val)) < lambda) {
            max_val <- max_val * 10
            if (max_val >= 1e9) stop(paste("\nERROR: Desired lambda is still",
                                           "not reached at max_val = 1e9"))
        }
    }

    if (direct_l) {
        op <- directL(op_fun, lower = 0, upper = max_val,
                      control = list(stopval = 0, xtol_rel = optim_tol))
        x_min <- op[["par"]]
        ldiff <- op[["value"]]
    } else {
        op <- optimize(op_fun, c(0, max_val), tol = optim_tol)
        x_min <- op[["minimum"]]
        ldiff <- op[["objective"]]
    }

    if (ldiff > precision)
        stop(paste0("\nERROR: Absolute difference between output lambda ",
                    "and desired lambda (", ldiff, ") exceeds required ",
                    "precision (", precision, ")."))

    L <- new_L_fun(x_min)

    return(L)

}





# Used in next two functions to make a list of arguments to pass
# to `adjust_lambda`:
make_al_arg_list <- function(L, lambda, adjust_part, al_args) {
    args <- list("L" = L, "lambda" = lambda, "adjust_part" = adjust_part)
    if (length(al_args) > 0) {
        if (is.null(names(al_args)) || any(names(al_args) == ""))
            stop("\nERROR: all items in al_args must be named")
        if (!all(names(al_args) %in% names(formals(adjust_lambda))))
            stop("\nERROR: all names in al_args must match args in adjust_lambda")
        for (n in names(args)) al_args[[n]] <- NULL # make sure these are ignored
        for (n in names(al_args)) args[[n]] <- al_args[[n]]
    }
    return(args)
}



#' Perturb Leslie matrix survivals stochastically
#'
#' @param L A numeric, square Leslie matrix.
#' @param sigma_s Standard deviation for the normal distribution used to
#'     generate deviates that are added to the logit-transformed survivals,
#'     before being back-transformed.
#' @param norm_rnd Single logical for whether random variables are normalized
#'     so that they always have a mean of zero. This results in
#'     survivals for the output Leslie matrix that always sum to the
#'     same number in logit-transformed space as those for the input matrix.
#'     This does NOT mean that the growth rate of the Leslie matrix will
#'     be the same.
#' @param norm_lambda Single logical for whether the output survivals
#'     are normalized so that the output Leslie matrix has the same growth rate
#'     as the input.
#'     Defaults to `FALSE`.
#' @param al_args Named list containing arguments to pass to [adjust_lambda()]
#'     that is used if `norm_lambda = TRUE`. The `adjust_lambda` arguments
#'     `L`, `lambda`, and `adjust_part` are ignored if provided here since
#'     they are defined inside this function.
#'     Defaults to `list()`.
#'
#' @returns A Leslie matrix with adjusted survivals.
#' @export
#'
random_survs <- function(L,
                         sigma_s,
                         norm_rnd = FALSE,
                         norm_lambda = FALSE,
                         al_args = list()) {
    if (sigma_s <= 0) return(L)
    lambda <- calc_lambda(L) # in case norm_lambda = TRUE
    survs <- L[row(L) - col(L) == 1]
    rnds <- rnorm(length(survs), 0, sigma_s)
    if (norm_rnd) rnds <- rnds - mean(rnds)
    survs <- inv_logit(logit(survs) + rnds)
    L[row(L) - col(L) == 1] <- survs
    if (norm_lambda) {
        args <- make_al_arg_list(L, lambda, "s", al_args)
        L <- do.call(adjust_lambda, args)
    }
    return(L)
}





#' Perturb Leslie matrix fecundities stochastically
#'
#' @param sigma_f Standard deviation for the normal distribution used to
#'     generate deviates that are exponentiated (i.e., `exp(x)` for random
#'     deviate(s) `x`), then multiplied by fecundities.
#' @param norm_rnd Single logical for whether random variables are
#'     normalized so that they always have a mean of one. This results in
#'     fecundities for the output Leslie matrix that always sum to the
#'     same number as those for the input matrix.
#'     This does NOT mean that the growth rate of the Leslie matrix will
#'     be the same.
#' @inheritParams random_survs
#'
#' @returns A Leslie matrix with adjusted fecundities.
#' @export
#'
random_fecunds <- function(L,
                           sigma_f,
                           norm_rnd = FALSE,
                           norm_lambda = FALSE,
                           al_args = list()) {
    if (sigma_f <= 0) return(L)
    lambda <- calc_lambda(L) # in case norm_lambda = TRUE
    rnds <- exp(rnorm(ncol(L), 0, sigma_f))
    if (norm_rnd) rnds <- rnds / mean(rnds)
    L[1,] <- L[1,] * rnds
    if (norm_lambda) {
        args <- make_al_arg_list(L, lambda, "f", al_args)
        L <- do.call(adjust_lambda, args)
    }
    return(L)
}
