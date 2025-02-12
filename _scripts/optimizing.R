
library(aphidsync)
library(gameofclones)




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
    op <- do.call(.optim, .args)
    return(op)
}



# ========================================================================*
# Main function ----
# ========================================================================*



winnowing_optim <- function(fn,
                            lower_bounds,
                            upper_bounds,
                            fn_args = list(),
                            box_control = list(maxit = 100, reltol = 1e-4),
                            fine_control = list(maxit = 500, reltol = 1e-6),
                            polished_control = list(maxit = 1000, reltol = 1e-8),
                            box_optimizer = optim,
                            fine_optimizer = optim,
                            polished_optimizer = optim,
                            n_bevals = 100L,
                            n_boxes = 1000L,
                            n_fine = 100L,
                            n_polished = 20L,
                            extra_optims = 0L) {


    stopifnot(is.function(fn))
    stopifnot(is.list(fn_args))
    stopifnot(length(lower_bounds) == length(upper_bounds))
    stopifnot(all(lower_bounds < upper_bounds))
    stopifnot(as.integer(n_bevals) == n_bevals)
    stopifnot(as.integer(n_boxes) == n_boxes)
    stopifnot(as.integer(n_fine) == n_fine)
    stopifnot(as.integer(n_polished) == n_polished)

    optim_pkgs <- list(box = packageName(environment(box_optimizer)),
                       fine = packageName(environment(fine_optimizer)),
                       polished = packageName(environment(polished_optimizer)))

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
        best_idx <- which(evals_i[,n_pars+1] == min(evals_i[,n_pars+1]))[[1]]
        best_pars <- evals_i[best_idx, 1:n_pars]
        op <- run_optim(box_optimizer, best_pars, fn, box_control, fn_args)

        box_evals[i,] <- c(best_pars, get_val(op, optim_pkgs$box))
    }

    # sort with lowest at the top:
    box_evals <- box_evals[order(box_evals[,n_pars+1]),]

    fines <- matrix(0.0, n_fine, ncol(box_evals))
    for (i in 1:n_fine) {
        pars <- box_evals[i,1:n_pars]
        op <- run_optim(fine_optimizer, pars, fn, fine_control, fn_args)
        fines[i,1:n_pars] <- pars
        fines[i,n_pars+1] <- get_val(op, optim_pkgs$fine)
    }
    # sort with lowest at the top:
    fines <- fines[order(fines[,n_pars+1]),]

    # keep all optimizer output in this case:
    polished_op <- lapply(1:n_polished, \(i) {
        pars <- fines[i,1:n_pars]
        op <- run_optim(polished_optimizer, pars, fn, polished_control, fn_args)
        return(op)
    })
    polished_vals <- sapply(polished_op, \(x) get_val(x, optim_pkgs$polished))
    best_op <- polished_op[[which(polished_vals == min(polished_vals))[[1]]]]
    not_converged <- !get_converged(best_op, optim_pkgs$polished)

    while (extra_optims > 0 && not_converged) {
        best_op <- run_optim(polished_optimizer, best_op$par, fn,
                             polished_control, fn_args)
        not_converged <- !get_converged(best_op, optim_pkgs$polished)
        extra_optims <- extra_optims - 1L
    }
    if (not_converged) warning("Final optimization has not converged.")

    return(best_op)

}






# ========================================================================*
# Test code ----
# ========================================================================*


.shape <- 2
.offset <- 0.1
.K <- 2000


# Susceptible line: no resistance, high population growth rate
line_s <- clonal_line("susceptible",
                      density_0 = 32,
                      surv_juv_apterous = "low",
                      surv_adult_apterous = "low",
                      repro_apterous = "low",
                      p_instar_smooth = 0)
# avoids weird rounding issue when `p_instar_smooth` = 0:
line_s$density_0 <- line_s$density_0 / sum(line_s$density_0) * 32
# Set starting abundances according to shape and offset above:
line_s$density_0[,1] <- beta_starts(.shape, .offset, sum(line_s$density_0),
                                   nrow(line_s$density_0))
line_s$density_0[,2] <- 0.0

sim_df <- gameofclones:::sim_gameofclones_full(line_s,
                                               n_fields = 1,
                                               n_plants = 1,
                                               max_t = 250,
                                               plant_check_gaps = 1,
                                               wilted_N = 1e6,
                                               mean_K = .K,
                                               clear_surv = 0,
                                               wasp_density_0 = 0,
                                               wasp_delay = 8,
                                               mum_smooth = 0,
                                               aphid_plant_disp_p = 0.1,
                                               plant_disp_mort = 0,
                                               wilted_mort = 0.3,
                                               extinct_N = 1,
                                               alate_b0 = -1e6,
                                               alate_b1 = 0) |>
    getElement("aphids") |>
    # group by aphid line and time (rep and field have only 1 option):
    filter(type != "mummy") |>
    group_by(time, line) |>
    summarize(N = sum(N), .groups = "drop")


re_df <- sim_df |>
    arrange(time) |>
    mutate(re = (lead(N) / N)^(1/(lead(time) - time)))



fn_args = list(L = line_s$leslie[,,1],
               re = re_df$re,
               time = re_df$time,
               max_shape = 1000)

# # Takes ~17 sec
# op <- winnowing_optim(fn = known_fit_aphids0,
#                       lower_bounds = c(  0, 0, 1e-9),
#                       upper_bounds = c(100, 1, 10e3),
#                       fn_args)
#
# op$par; c(.shape, .offset, .K)
# op$par - c(.shape, .offset, .K)


cntrl <- list(box = list(maxit = 100, reltol = 1e-4),
              fine = list(maxit = 500, reltol = 1e-6),
              polished = list(maxit = 1000, reltol = 1e-8))
# for nloptr
cntrl2 <- list(box = list(maxeval = 100, xtol_rel = 1e-4),
              fine = list(maxeval = 500, xtol_rel = 1e-6),
              polished = list(maxeval = 1000, xtol_rel = 1e-8))
# for minqa
cntrl3 <- list(box = list(maxfun = 100),
              fine = list(maxfun = 500),
              polished = list(maxfun = 1000))



foo1 <- function(.t) do.call(optim, c(list(par = c(50, 0.5, 1000),
                                           fn = known_fit_aphids0,
                                           control = cntrl[[.t]]),
                                      fn_args))
foo2 <- function(.t) do.call(nloptr::bobyqa, c(list(x0 = c(50, 0.5, 1000),
                                                    fn = known_fit_aphids0,
                                                    control = cntrl2[[.t]]),
                                               fn_args))
foo3 <- function(.t) do.call(nloptr::newuoa, c(list(x0 = c(50, 0.5, 1000),
                                                    fn = known_fit_aphids0,
                                                    control = cntrl2[[.t]]),
                                               fn_args))
foo4 <- function(.t) do.call(minqa::bobyqa, c(list(par = c(50, 0.5, 1000),
                                                   fn = known_fit_aphids0,
                                                   control = cntrl3[[.t]]),
                                              fn_args))
foo5 <- function(.t) do.call(minqa::newuoa, c(list(par = c(50, 0.5, 1000),
                                                   fn = known_fit_aphids0,
                                                   control = cntrl3[[.t]]),
                                              fn_args))



microbenchmark::microbenchmark(optim = foo1("box"), nl_bobyqa = foo2("box"),
                               nl_newuoa = foo3("box"), mi_bobyqa = foo4("box"),
                               mi_newuoa = foo5("box"))
