

# ========================================================================*
# Helper functions ----
# ========================================================================*

# Get the value of the objective at the best set of parameters found.
get_val <- function(.op, .pkg) {
    if (.pkg %in% c("stats", "nloptr")) {
        val <- .op$value
    } else {
        val <- .op$fval
    }
    return(val)
}
# Get whether an optimization has converged.
get_converged <- function(.op, .pkg) {
    if (.pkg == "minqa") {
        conv <- .op$ierr == 0L
    } else if (.pkg == "nloptr") {
        conv <- .op$convergence > 0L
    } else {
        conv <- .op$convergence == 0L
    }
    return(conv)
}
# Combine arguments and run optimizer
run_optim <- function(.optim, .pars, .fn, .control, .fn_args) {
    .args <- c(list(par = .pars, fn = .fn, control = .control), .fn_args)
    # Adjust arguments across packages:
    if (packageName(environment(.optim)) == "nloptr") {
        if (!is.null(.args[["control"]][["maxit"]])) {
            .args[["control"]][["maxeval"]] <- .args[["control"]][["maxit"]]
            .args[["control"]][["maxit"]] <- NULL
        }
        if (!is.null(.args[["control"]][["reltol"]])) {
            .args[["control"]][["xtol_rel"]] <- .args[["control"]][["reltol"]]
            .args[["control"]][["reltol"]] <- NULL
        }
        .args[["x0"]] <- .pars
        .args[["par"]] <- NULL
    }
    if (packageName(environment(.optim)) == "minqa") {
        if (!is.null(.args[["control"]][["maxit"]])) {
            .args[["control"]][["maxfun"]] <- .args[["control"]][["maxit"]]
            .args[["control"]][["maxit"]] <- NULL
        }
        # No obvious equivalent in `minqa`, so just remove this:
        .args[["control"]][["reltol"]] <- NULL
    }
    # `lower` and `upper` never belong in the control argument, but sometimes
    # they can be used in the optimization function itself (depending on the
    # algorithm):
    for (n in c("lower", "upper")) {
        if (!is.null(.args[["control"]][[n]])) {
            if (n %in% names(formals(.optim))) {
                .args[[n]] <- .args[["control"]][[n]]
            }
            .args[["control"]][[n]] <- NULL
        }
    }
    op <- do.call(.optim, .args)
    return(op)
}



# ========================================================================*
# Main function ----
# ========================================================================*





#' Winnowing optimization
#'
#' @param fn Objective function to minimize.
#' @param lower_bounds Lower bounds of boxes for each parameter.
#' @param upper_bounds Upper bounds of boxes for each parameter.
#' @param fn_args List containing other arguments to use for `fn`.
#'     Defaults to `list()`.
#' @param box_control List containing arguments to use for the `control`
#'     argument for the optimization that occurs for each box.
#'     Defaults to `list(maxit = 100, reltol = 1e-4)`
#'     (these names are adjusted if using an optimizer other than `optim`).
#' @param fine_control List containing arguments to use for the `control`
#'     argument for the fine optimizations.
#'     Defaults to `list(maxit = 500, reltol = 1e-6)`
#'     (these names are adjusted if using an optimizer other than `optim`).
#' @param polished_control List containing arguments to use for the `control`
#'     argument for the polished optimizations.
#'     Defaults to `list(maxit = 1000, reltol = 1e-8)`
#'     (these names are adjusted if using an optimizer other than `optim`).
#' @param box_optim Function to use for the optimization for each box.
#'     Must be `stats::optim` or a function from packages `minqa` or `nloptr`.
#'     Defaults to `stats::optim`.
#' @param fine_optim Function to use for the fine optimizations.
#'     Must be `stats::optim` or a function from packages `minqa` or `nloptr`.
#'     Defaults to `stats::optim`.
#' @param polished_optim Function to use for the polished optimizations.
#'     Must be `stats::optim` or a function from packages `minqa` or `nloptr`.
#'     Defaults to `stats::optim`.
#' @param n_bevals Number of evaluations of `fn` per box. Defaults to `100L`.
#' @param n_boxes Number of boxes. Defaults to `1000L`.
#' @param n_fine Number of fine optimizations.
#'     Must be `> n_boxes` and `< n_polished`.
#'     Defaults to `100L`.
#' @param n_polished Number of polished optimizations.
#'     Must be `> n_fine`.
#'     Defaults to `20L`.
#' @param n_outputs Number of top output object(s) to return.
#'     Must be `> n_polished`. Defaults to `3L`.
#' @param extra_optims Number of potential extra rounds of optimizations to run
#'     if the best polished optimization still hasn't converged.
#'     Defaults to `10L`.
#'
#' @importFrom stats optim
#'
#' @returns `n_outputs` object(s) of the class returned by `polished_optim`.
#'
#'
#' @export
#'
winnowing_optim <- function(fn,
                            lower_bounds,
                            upper_bounds,
                            fn_args = list(),
                            box_control = list(maxit = 100, reltol = 1e-4),
                            fine_control = list(maxit = 500, reltol = 1e-6),
                            polished_control = list(maxit = 1000, reltol = 1e-8),
                            box_optim = optim,
                            fine_optim = optim,
                            polished_optim = optim,
                            n_bevals = 100L,
                            n_boxes = 1000L,
                            n_fine = 100L,
                            n_polished = 20L,
                            n_outputs = 3L,
                            extra_optims = 10L) {


    # Type and length checks:
    stopifnot(is.function(fn))
    stopifnot(is.numeric(lower_bounds) && is.numeric(upper_bounds))
    stopifnot(length(lower_bounds) == length(upper_bounds))
    stopifnot(is.list(fn_args))
    stopifnot(is.list(box_control))
    stopifnot(is.list(fine_control))
    stopifnot(is.list(polished_control))
    stopifnot(is.function(box_optim))
    stopifnot(is.function(fine_optim))
    stopifnot(is.function(polished_optim))
    stopifnot(length(n_bevals) == 1L && as.integer(n_bevals) == n_bevals)
    stopifnot(length(n_boxes) == 1L && as.integer(n_boxes) == n_boxes)
    stopifnot(length(n_fine) == 1L && as.integer(n_fine) == n_fine)
    stopifnot(length(n_polished) == 1L && as.integer(n_polished) == n_polished)
    stopifnot(length(n_outputs) == 1L && as.integer(n_outputs) == n_outputs)
    stopifnot(length(extra_optims) == 1L && as.integer(extra_optims) == extra_optims)

    # Value checks:
    stopifnot(all(lower_bounds < upper_bounds))
    stopifnot(n_bevals > 0L)
    stopifnot(n_boxes > 0L)
    stopifnot(n_fine > 0L && n_fine < n_boxes)
    stopifnot(n_polished > 0L && n_polished < n_fine)
    stopifnot(n_outputs > 0L && n_outputs < n_polished)
    stopifnot(extra_optims >= 0L)


    optim_pkgs <- list(box = packageName(environment(box_optim)),
                       fine = packageName(environment(fine_optim)),
                       polished = packageName(environment(polished_optim)))

    if (!all(optim_pkgs %in% c("stats", "minqa", "nloptr"))) {
        stop(paste("this function only programmed for stats::optim and",
                   "optimizers from minqa and nloptr packages"))
    }

    mids <- (upper_bounds + lower_bounds) / 2
    steps <- (upper_bounds - mids) / n_boxes
    n_pars <- length(mids)

    box_evals <- matrix(0.0, n_boxes, n_pars+1)
    evals_i <- matrix(0.0, n_bevals, n_pars+1)

    for (i in 1:n_boxes) {
        for (j in 1:n_bevals) {
            .pars <- runif(n_pars, mids - steps * i, mids + steps * i)
            .val <- do.call(fn, c(list(.pars), fn_args))
            evals_i[j,] <- c(.pars, .val)
        }
        min_evals <- min(evals_i[,n_pars+1], na.rm = TRUE)
        best_idx <- which(evals_i[,n_pars+1] == min_evals)[[1]]
        best_pars <- evals_i[best_idx, 1:n_pars]
        op <- run_optim(box_optim, best_pars, fn, box_control, fn_args)

        box_evals[i,] <- c(best_pars, get_val(op, optim_pkgs$box))
    }

    # sort with lowest at the top:
    box_evals <- box_evals[order(box_evals[,n_pars+1]),]

    fines <- matrix(0.0, n_fine, ncol(box_evals))
    for (i in 1:n_fine) {
        pars <- box_evals[i,1:n_pars]
        op <- run_optim(fine_optim, pars, fn, fine_control, fn_args)
        fines[i,1:n_pars] <- pars
        fines[i,n_pars+1] <- get_val(op, optim_pkgs$fine)
    }
    # sort with lowest at the top:
    fines <- fines[order(fines[,n_pars+1]),]

    # keep all optimizer output in this case:
    polished_op <- lapply(1:n_polished, \(i) {
        pars <- fines[i,1:n_pars]
        op <- run_optim(polished_optim, pars, fn, polished_control, fn_args)
        return(op)
    })
    polished_vals <- sapply(polished_op, \(x) get_val(x, optim_pkgs$polished))
    sort_idx <- order(polished_vals)
    polished_op <- polished_op[sort_idx]
    polished_vals <- polished_vals[sort_idx]

    best_ops <- polished_op[1:n_outputs]
    not_conv <- sapply(best_ops, \(o) !get_converged(o, optim_pkgs$polished))

    for (i in 1:n_outputs) {
        eo <- as.integer(extra_optims)
        while (eo > 0L && not_conv[i]) {
            best_ops[[i]] <- run_optim(polished_optim, best_ops[[i]]$par, fn,
                                       polished_control, fn_args)
            not_conv[i] <- !get_converged(best_ops[[i]], optim_pkgs$polished)
            eo <- eo - 1L
        }
        if (not_conv[i]) warning("Final optimization ", i, " has not converged.")
    }

    return(best_ops)

}



