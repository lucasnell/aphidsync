

# Function to convert from scale of interest to that for fitting:
trans_pars <- function(pars, trans_base) {
    if (!length(pars) %in% 2:4) stop("length(pars) should be 2, 3, or 4")
    pars[[1L]] <- log(pars[[1L]] - 1.0, trans_base)
    pars[[2L]] <- logit(pars[[2L]])
    if (length(pars) > 2) pars[[3L]] <- log(pars[[3L]], trans_base)
    return(pars)
}

# Function to convert from scale of fitting to that of interest:
# (I'm not transforming 4th parameter (if there is one) bc its use isn't
# straightforward; see description above `known_fit_aphids0` for more info.)
inv_trans_pars <- function(pars, trans_base) {
    if (!length(pars) %in% 2:4) stop("length(pars) should be 2, 3, or 4")
    pars[[1L]] <-trans_base^pars[[1L]] + 1.0
    pars[[2L]] <- inv_logit(pars[[2L]])
    if (length(pars) > 2) pars[[3L]] <-trans_base^pars[[3L]]
    return(pars)
}



#' Fit starting abundance for known Leslie matrix (optionally with uncertainty)
#'
#' @param sim_df A dataframe containing the time and abundance columns.
#'     Requires columns `"time"` and `"N"`.
#' @param L Leslie matrix used for fitting.
#' @param N0 Starting abundance.
#' @param K Density dependence if known, or `NA` if unknown. Defaults to `NA`.
#' @param adjust_L Single integer for whether and how to fit parameters
#'     that adjust the Leslie matrix.
#'     A value of `0L` results in no Leslie adjustment,
#'     `1L` adjusts the Leslie using 1 parameter (i.e., `L[,1] = L[,1] * fecund_x`),
#'     and `2L` adjusts the Leslie using separate parameters for
#'     the fecundities and survivals.
#'     Fecundities (`f`) are adjusted as `f = f * fecund_x`.
#'     Survivals (`s`) are adjusted via `s = inv_logit(logit(s) + surv_x)`.
#'     Defaults to `0L`.
#' @param max_shape Single numeric for the maximum value for the shape
#'     parameter. Depending on the number of stages, greater increases to
#'     the shape parameter have little effect. Defaults to `800`.
#' @param compare_N Single logical for whether to use `N` or per-capita growth
#'     when comparing between input data and object function simulations.
#' @param fit_max_t Single integer for the number of time steps to use
#'     for estimating the agreement between observed vs fitted values.
#'     Too high of a value here results in density dependence coming into
#'     play more strongly. Defaults to `30L`.
#' @param est_K_start Single integer for the time at which to start estimating
#'     density dependence. Defaults to `50L`.
#' @param est_K_min_n Single integer for the number of time steps to use
#'     for estimating density dependence. Defaults to `50L`.
#' @param max_fecund_x Single numeric for the maximum value of `fecund_x` to
#'     search for within `winnowing_optim`.
#'     Parameter `fecund_x` is the value that gets multiplied by
#'     fecundities when adjusting the Leslie matrix using
#'     2 parameters (i.e., `adjust_L = 2L`). Ignored if `adjust_L != 2L`.
#'     Defaults to `10`.
#' @param max_surv_x Single numeric for the maximum absolute value of `surv_x`
#'     to search for within `winnowing_optim`.
#'     Parameter `surv_x` is the value that gets added to
#'     logit-transformed survivals when adjusting the Leslie matrix using
#'     2 parameters (i.e., `adjust_L = 2L`). Ignored if `adjust_L != 2L`.
#'     Defaults to `10`.
#' @param trans_base Single numeric indicating the base to use for
#'     transforming the shape and (optionally) `fecund_x` parameters when
#'     fitting.
#'     Each parameter `z` will be transformed as `trans_base^z` inside the
#'     objective function, and back-transformed as `log(z, trans_base)`
#'     in this one.
#'     Using `exp(1)` here seems to cause the optimizer to have a hard time.
#'     The default for this is `1.36`, which is about `exp(1) / 2`.
#' @param optim_args List containing arguments to pass to `winnowing_optim`.
#'     Defaults to `list()`.
#' @param return_optims Single logical for whether to return all polished
#'     optimization objects from `winnowing_optim`.
#'     The alternative is to just return the fitted parameter values with
#'     a vector indicating the greatest divergence between any of the polished
#'     fits for each parameter.
#'     Defaults to `FALSE`.
#'
#' @returns A named vector indicating the parameter values from the fit.
#'
#' @export
#'
fit_known_leslie <- function(sim_df, L, N0,
                             K = NA_real_,
                             adjust_L = 0L,
                             max_shape = 800,
                             compare_N = FALSE,
                             fit_max_t = 30L,
                             est_K_start = 50L,
                             est_K_min_n = 50L,
                             max_fecund_x = 10,
                             max_surv_x = 10,
                             trans_base = 1.36,
                             optim_args = list(),
                             return_optims = FALSE) {

    # Used for making numbers nearly zero but not quite:
    VERY_SMALL <- max(.Machine$double.eps, .Machine$double.neg.eps)

    type_checker(sim_df, "sim_df", "data.frame", l = NA)
    stopifnot(all(c("time", "N") %in% colnames(sim_df)))
    type_checker(L, "L", "matrix", l = NA, .min = 0)
    stopifnot(nrow(L) == ncol(L))
    type_checker(N0, "N0", "numeric", .min = VERY_SMALL)
    if (!isTRUE(is.na(K))) type_checker(K, "K", "numeric", .min = VERY_SMALL)
    type_checker(adjust_L, "adjust_L", "integer", .min = 0L, .max = 2L)
    type_checker(max_shape, "max_shape", "numeric", .min = 1)
    type_checker(compare_N, "compare_N", "logical")
    type_checker(fit_max_t, "fit_max_t", "integer")
    type_checker(est_K_start, "est_K_start", "integer", .min = 1L)
    type_checker(est_K_min_n, "est_K_min_n", "integer", .min = 1L)
    type_checker(max_fecund_x, "max_fecund_x", "numeric", .min = 2 * VERY_SMALL)
    # note on line below: 1+VERY_SMALL is too small!
    type_checker(trans_base, "trans_base", "numeric", .min = 1 + 1e-6)
    type_checker(optim_args, "optim_args", "list", l = NA)
    type_checker(return_optims, "return_optims", "logical")

    if (nrow(sim_df) < fit_max_t) {
        stop(sprintf("sim_df has %i rows and fit_max_t = %i",
                     nrow(sim_df), fit_max_t))
    }
    if (is.na(K) && sum(sim_df$time > est_K_start) < est_K_min_n) {
        stop(sprintf(paste("sim_df has %i rows with time > est_K_start",
                           "and needs %i to estimate K, so get more data or",
                           "change est_K_start or est_K_min_n."),
                     sum(sim_df$time > est_K_start), est_K_min_n))
    }
    if ("line" %in% colnames(sim_df) && length(unique(sim_df$line)) > 1L) {
        stop("sim_df cannot have >1 line")
    }

    time <- sim_df$time[1:fit_max_t]
    stopifnot(all(time >= 0))
    if (compare_N) {
        observed <- sim_df$N[1:fit_max_t]
    } else {
        re_df <- sim_df |>
            arrange(time) |>
            mutate(re = (lead(N) / N)^(1/(lead(time) - time)))
        observed <- re_df$re[1:fit_max_t]
    }

    # If not provided as a known value, estimate K from time series and
    # Leslie matrix:
    if (is.na(K)) {
        # Note: simulations showed that `mean(N)` (versus `median(N)`,
        # `exp(mean(log(N)))`, or `exp(median(log(N)))`) works best for
        # this prediction. This is true also with measurement error.
        K <- mean(sim_df$N[sim_df$time > est_K_start]) /
            (max(abs(eigen(L)[["values"]])) - 1)
        if (is.na(K)) stop("is.na(K) after prediction attempt")
    }

    # pointer to a c++ object used in objective function:
    fit_ptr <- make_known_fit_ptr(K = K, L = L, obs = observed, time = time,
                                  max_shape = max_shape, N0 = N0,
                                  compare_N = compare_N, trans_base = trans_base)


    # Do optimization using winnowing approach
    op_args <- optim_args
    op_args[["fn"]] <- known_fit_aphids0
    op_args[["lower_bounds"]] <- c(rep(VERY_SMALL, 3) + c(1,0,0), -max_surv_x) |>
        base::`[`(1:(2L+adjust_L)) |>
        trans_pars(trans_base = trans_base)
    op_args[["upper_bounds"]] <- c(max_shape, 1-VERY_SMALL,
                                   max_fecund_x, max_surv_x) |>
        base::`[`(1:(2L+adjust_L)) |>
        trans_pars(trans_base = trans_base)
    op_args[["fn_args"]] <- list(fit_info_ptr = fit_ptr)


    ops <- do.call(winnowing_optim, op_args)
    # Reverse transformations that occur inside `known_fit_aphids0`
    # (I'm not transforming 4th parameter if Leslie fitting occurred
    #  bc its use isn't straightforward; see description above `known_fit_aphids0`
    #  for more info):
    for (i in 1:length(ops)) {
        if (!"par" %in% names(ops[[i]])) stop("unknown output parameter name")
        ops[[i]][["par"]] <- inv_trans_pars(ops[[i]][["par"]], trans_base)
    }

    if (return_optims) return(ops)

    op <- ops[[1]]

    pars <- c(op$par, K, width99(op$par[[1]]),
              med_age(op$par[[1]], op$par[[2]], ncol(L)))
    if (adjust_L == 0L) {
        names(pars) <- c("shape", "offset", "K", "width99", "med_age")
    } else if (adjust_L == 1L) {
        names(pars) <- c("shape", "offset", "fecund_x", "K", "width99", "med_age")
    } else {
        names(pars) <- c("shape", "offset", "fecund_x", "surv_x", "K", "width99", "med_age")
    }

    # This will inform whether the different optimizations returned from
    # `winnowing_optim` are the same (or roughly so)
    if (length(ops) > 1) {
        par_diffs <- sapply(ops[-1], \(o) {
            sh <- o$par[[1]]
            of <- o$par[[2]]
            o_pars <- c(o$par, K, width99(sh), med_age(sh, of, ncol(L)))
            names(o_pars) <- names(pars)
            return(o_pars - pars)
        }) |>
            apply(1, \(x) { x[which(abs(x) == max(abs(x)))[1]] })
        if (any(is.na(par_diffs))) stop("any(is.na(par_diffs))")
    } else par_diffs <- pars * 0 # this way keeps names
    attr(pars, "par_diffs") <- par_diffs

    return(pars)

}








