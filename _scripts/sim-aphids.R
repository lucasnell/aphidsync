
suppressPackageStartupMessages({
    library(tidyverse)
    library(viridisLite)
    library(gameofclones)
    library(aphidsync)
    library(future.apply)
    library(progressr)
    library(patchwork)
})

plan(multisession, workers = max(parallel::detectCores() - 2L, 1L))
handlers(global = TRUE)
handlers("progress")


# ggplot2 theme:
theme_set(theme_classic() %+replace%
              theme(strip.background = element_blank(),
                    strip.text = element_text(size = 11),
                    legend.background = element_blank(),
                    plot.title = element_text(size = 14, hjust = 0)))

if (interactive() && file.exists(".Rprofile")) source(".Rprofile")



col_pal <- c(resistant = viridis(100)[50],
             susceptible = viridis(100)[95],
             wasp = viridis(100)[1])

# Define clonal lines. Don't re-define these!
# Susceptible line: no resistance, high population growth rate
line_s <- clonal_line("susceptible",
                      density_0 = 32,
                      # surv_juv_apterous = "high",
                      # surv_adult_apterous = "high",
                      # repro_apterous = "high")
                      surv_juv_apterous = "low",
                      surv_adult_apterous = "low",
                      repro_apterous = "low",
                      p_instar_smooth = 0)
# avoids weird rounding issue when `p_instar_smooth` = 0:
line_s$density_0 <- line_s$density_0 / sum(line_s$density_0) * 32

# # Resistant line: high resistance, low population growth rate
# line_r <- clonal_line("resistant",
#                       density_0 = 32,
#                       resistant = TRUE,
#                       surv_paras = 0.57,
#                       surv_juv_apterous = "low",
#                       surv_adult_apterous = "low",
#                       repro_apterous = "low",
#                       p_instar_smooth = 0)
# for(i in 1:3) line_r$leslie[1,,i] <- 0.8 * line_r$leslie[1,,i]



beta_startsR <- function(shape, offset, total_aphids0, compartments = 29L) {

    stopifnot(shape > 0)
    stopifnot(offset >= 0 && offset <= 1)
    stopifnot(compartments > 0 && total_aphids0 > 0)
    stopifnot(compartments == round(compartments))
    stopifnot(total_aphids0 == round(total_aphids0))
    times <- seq(0, 1, length.out = round(compartments + 1)) + offset
    correct <- \(x, f) ifelse(times > 1, f(x, 1), x)
    pbeta_vals <- pbeta(correct(times, `-`), shape1 = shape, shape2 = shape) |>
        correct(`+`)
    aphids0 <- total_aphids0 * diff(pbeta_vals)
    if(round(sum(aphids0)) != round(total_aphids0)){
        warning(paste("beta_startsR magnitude error (",
                      round(sum(aphids0)) - round(total_aphids0), ", shape = ",
                      round(shape, digits = 2), ", offset = ",
                      round(offset, digits = 2), sep = ''), immediate. = TRUE)
    }
    if(round(length(aphids0)) != round(compartments)){
        warning(paste("beta_startsR length error, shape = ",
                      round(shape, digits = 2), ", offset = ",
                      offset, sep=''), immediate. = TRUE)
    }
    return(aphids0)
}






sim_aphids <- function(.shape, .offset, .K, sample_filter = FALSE) {

    line <- line_s
    line$density_0[,1] <- beta_starts(.shape, .offset, sum(line$density_0),
                                      nrow(line$density_0))
    line$density_0[,2] <- 0.0

    sim_df <- gameofclones:::sim_gameofclones_full(line,
                                                   n_fields = 1,
                                                   # n_plants = 16,
                                                   n_plants = 1,
                                                   max_t = 250,
                                                   # aphid_demog_error = TRUE,
                                                   # environ_error = TRUE,
                                                   # >>>>>>>>>>>>>>
                                                   # plant_check_gaps = c(3, 4),
                                                   # extra_plant_removals =
                                                   #     cbind(8, 4),
                                                   # wilted_N = 400,
                                                   # >>>>>>>>>>>>>>
                                                   plant_check_gaps = 1,
                                                   wilted_N = 1e6,
                                                   # >>>>>>>>>>>>>
                                                   mean_K = .K,  # 1800,
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

    if (sample_filter){
        sim_df <- sim_df |>
            # Filter for sampling twice per week.
            # Starting with the first day, sample that day of the week and
            # one other three days later:
            filter(time %% 7L == 0 | time %% 7L == 3L)
    }

    return(sim_df)

}







#
# Use within optim to find shape parameters that minimize the differences
# between observed and predicted Re values.
#
fit_aphids0R <- function(pars, known_L, re_df, fecund, match_lambda, fit_survs,
                         .plot = FALSE, L_upper_bound = 1000, max_shape = 800) {

    # pars <- c(533.764, 393.973)

    .K <- 1800

    if (any(pars < 0) || pars[2] > 1) return(1e10)
    .shape <- pars[[1]]
    .offset <- pars[[2]]
    # .K <- pars[[3]]
    if (.shape > max_shape) return(1e10)

    aphids0 <- beta_starts(.shape, .offset, total_aphids0 = 32)

    if (known_L) {
        L <- line_s$leslie[,,1]
    } else {
        stopifnot(length(pars) >= 4)
        if (any(pars[3:4] == 0)) return(1e10)
        if (fit_survs && length(pars) < 5) stop("fit_survs && length(pars) < 5")
        if (pars[5] > 1) return(1e10)
        n_stages <- 29L
        adult_stage <- 9L # first adult stage
        # Setup starting leslie matrix with only survivals = 1 and with
        # fecundities summing to 1. Lambda is also 1.
        L <- make_L1(pars[[3]], pars[[4]])
        if (fit_survs) {
            stages0 <- -(n_stages-2L):0
            L[row(L) - col(L) == 1] <- 1 - pars[[5]] / sqrt(1 + stages0^2)
        }

        if (match_lambda) {

            lambda <- fecund

            # If inputs exceed bounds of regression used to fit L, then fit
            # using an optimizer (optimizing is ~500x slower)
            if (pars[[3]] > 5 || pars[[4]] > 20 || lambda > 1.6 || lambda < 1.05) {
                L_op <- optimize(\(x) {
                    Lx <- L
                    Lx[1, 9:29] <- Lx[1, 9:29] * x
                    lambda1 <- max(abs(eigen(Lx, only.values = TRUE)[["values"]]))
                    return(abs(lambda1 - lambda))
                }, c(0, L_upper_bound))
                x <- L_op$minimum
            } else {
                x <- calc_L(pars[[3]], pars[[4]], lambda)[1,1]
            }

        } else {

            x <- fecund / max(L[1, adult_stage:n_stages])

        }

        L[1, adult_stage:n_stages] <- L[1, adult_stage:n_stages] * x

    }

    # Extrapolate through time:
    Nmat  <- matrix(NA_real_, nrow = max(re_df$time) + 1, ncol = length(aphids0))
    Nmat[1, ] <- aphids0
    for (i in 2:(nrow(Nmat))) {
        Ni <- t(Nmat[(i-1),,drop=FALSE])
        z <- sum(Ni)
        S <- 1 / (1 + z / .K)
        Nmat[i, ] <- S * t(L %*% Ni)
    }

    sumNvec <- apply(Nmat, 1, sum)
    re_pred <- (lead(sumNvec) / sumNvec)[(re_df$time + 1L)]

    if (.plot) {
        plot(re_pred, ylim = c(0, max(c(re_pred, re_df$re), na.rm = TRUE)),
             type = "l", ylab = "Re predicted (black) and observed (red)",
             xlab = "time", lwd = 2,
             main = sprintf("shape = %.4f, offset = %.4f", .shape, .offset))
        lines(re_df$re, col = "red", lwd = 2)
    }

    sae <- sum(abs(re_pred - re_df$re), na.rm = TRUE)
    if (is.na(sae) | sae > 1e10) sae <- 1e10

    return(sae)

}







#
# Fit the two shape parameters for a set of simulations.
# NOTE: This version doesn't allow for multiple lines (yet).
#
fit_sims <- function(sim_df, cycles = TRUE, known_L = TRUE, cpp = TRUE,
                     match_lambda = FALSE, fit_survs = TRUE,
                     fit_L_control = list(), method = "optim") {

    stopifnot(length(unique(sim_df$line)) == 1L)

    method <- match.arg(method, c("optim", "BOBYQA", "COBYLA", "NELDERMEAD", "SBPLX"))

    if (length(fit_L_control) > 0) {
        stopifnot(is.list(fit_L_control))
        stopifnot(!is.null(names(fit_L_control)) && all(names(fit_L_control) != ""))
    }

    if (nrow(sim_df) < 3L) return("low rows") # special meaning in `test_sim`

    # n_stages <- formals(beta_starts)[["compartments"]]
    # stopifnot(n_stages %in% c(14, 29))
    # adult_stage <- ifelse(n_stages == 14, 5, 9) # first adult stage
    #
    # # Setup starting leslie matrix with only survivals = 1 (without fecundities)
    # leslie <- matrix(0, n_stages, n_stages)
    # # survivals (set to 1 for now):
    # leslie[(row(leslie) - col(leslie) == 1)] <- 1
    # # fecundities (to be fit):
    # leslie[1, adult_stage:n_stages] <- 1
    # # If you want to test known fecundities:
    # # leslie[1, adult_stage:n_stages] <- line_s$leslie[1, adult_stage:n_stages,1]

    # `fecund` parameter depends on whether we're trying to match Leslie matrix
    # lambda or max fecundity
    lambda <- max(abs(eigen(line_s$leslie[,,1], only.values = TRUE)[["values"]]))
    if (match_lambda) {
        fecund <- lambda
    } else fecund <- max(line_s$leslie[,,1])


    re_df <- sim_df |>
        arrange(time) |>
        mutate(re = (lead(N) / N)^(1/(lead(time) - time)))

    # amr_df <- sim_df
    # if (cycles) {
    #     amr_df <- amr_df |>
    #         filter(time %% 29L == 0)
    # }
    # amr_df <- amr_df |>
    #     mutate(amr = (lead(N) / N)^(1/(lead(time) - time))) |>
    #     filter(!is.na(amr))
    #
    #
    # amr_f <- approxfun(amr_df$time, y = amr_df$amr, yleft = amr_df$amr[1],
    #                    yright = tail(amr_df$amr, 1))
    #
    # amr_pred = amr_f(1:max(re_df$time))
    #
    # if (known_L) {
    #     L <- line_s$leslie[,,1]
    # } else {
    #     L_fit_fun <- function(x, .amr) {
    #         L <- leslie
    #         L[1, adult_stage:n_stages] <- x[1]
    #         lam <- max(abs(eigen(L, only.values = TRUE)[["values"]]))
    #         return(abs(lam - .amr))
    #     }
    #     L_op <- optim(1, L_fit_fun, .amr = amr)
    #
    #     Lvec <- lapply(amr_pred, \(amr) {
    #         if (op$convergence != 0) stop("Leslie could not converge")
    #         L <- leslie
    #         L[1, adult_stage:n_stages] <- op$par[1]
    #         return(L)
    #     })
    # }

    # Fit across multiple starting offsets, then choose best-fitting one:
    val <- 1e10
    op <- NULL
    regr_ptr <- NULL
    fit_L_args <- list(L_method = "BOBYQA",
                     L_upper_bound = 1000,
                     L_tol = 1e-4,
                     L_max_iters = 1000,
                     max_shape = 800)
    optimize_f <- NULL
    if (method == "optim") {
        optimize_f <- optim
    } else if (method == "BOBYQA") {
        optimize_f <- nloptr::bobyqa
    } else if (method == "COBYLA") {
        optimize_f <- nloptr::cobyla
    } else if (method == "NELDERMEAD") {
        optimize_f <- nloptr::neldermead
    } else if (method == "SBPLX") {
        optimize_f <- nloptr::sbplx
    } else stop("unknown method")
    if (cpp) {
        regr_ptr <- make_regr_ptr()
    }
    if (length(fit_L_control) > 0) {
        all(names(fit_L_control) %in% names(fit_L_args))
        for (n in names(fit_L_control)) fit_L_args[[n]] <- fit_L_control[[n]]
    }
    n_starts <- 100L
    if (known_L) {
        known_L_mat <- line_s$leslie[,,1]
        guess_list <- lapply(1:n_starts, \(i) runif(2, 0, c(10, 1)))
    } else {
        known_L_mat <- matrix(0, 0, 0)
        if (fit_survs) {
            guess_list <- lapply(1:n_starts, \(i) runif(5, 0, c(10, 1, 10, 20, 1)))
        } else {
            guess_list <- lapply(1:n_starts, \(i) runif(4, 0, c(10, 1, 10, 20)))
        }
    }
    # Choose 10 random reps to have very low and high offset values.
    # This helps with fitting.
    n_hilo <- 10L
    off_inds <- sample.int(n_starts, n_hilo)
    for (i in off_inds[1:(n_hilo %/% 2L)])
        guess_list[[i]][2] <- runif(1, 0, 1e-6)
    for (i in off_inds[(n_hilo %/% 2L + 1L):n_hilo])
        guess_list[[i]][2] <- runif(1, 1 - 1e-6, 1)
    # Similar for survival parameter (if used):
    if (!known_L && fit_survs) {
        spar_inds <- sample.int(n_starts, n_hilo)
        for (i in spar_inds[1:(n_hilo %/% 2L)])
            guess_list[[i]][5] <- runif(1, 0, 1e-6)
        for (i in spar_inds[(n_hilo %/% 2L + 1L):n_hilo])
            guess_list[[i]][5] <- runif(1, 1 - 1e-6, 1)
    }
    for (s in 1:n_starts) {
        guess <- guess_list[[s]]
        if (cpp) {
            new_op <- optimize_f(guess, fit_aphids0,
                            known_L = known_L_mat,
                            re = re_df$re, time = re_df$time,
                            fecund = fecund, regr_ptr = regr_ptr,
                            match_lambda = match_lambda,
                            fit_survs = fit_survs,
                            max_shape = fit_L_args[["max_shape"]],
                            L_method = fit_L_args[["L_method"]],
                            L_upper_bound = fit_L_args[["L_upper_bound"]],
                            L_tol = fit_L_args[["L_tol"]],
                            L_max_iters = fit_L_args[["L_max_iters"]])
        } else {
            new_op <- optimize_f(guess, fit_aphids0R, fecund = fecund,
                            match_lambda = match_lambda,
                            fit_survs = fit_survs,
                            known_L = known_L, re_df = re_df,
                            L_upper_bound = fit_L_args[["L_upper_bound"]],
                            max_shape = fit_L_args[["max_shape"]])
        }
        if (new_op$convergence == 0 && new_op$value < val) {
            val <- new_op$value
            op <- new_op
        }
        if (new_op$convergence == 0 && !is.null(op) && new_op$value == val) {
            if (!isTRUE(all.equal(op$par, new_op$par))) {
                # This triggers a specific error in `test_sim`:
                return("multiple models")
            }
        }
    }


    if (is.null(op)) return("no converged") # specific meaning in `test_sim`

    pars <- op$par
    names(pars) <- c("shape", "offset", "wshape", "wscale", "spar")[1:length(op$par)]

    return(pars)

}



# Width of 99th quartile (offset ignored):
width99 <- function(shape, offset) {
    sapply(shape, \(sh) {
        q99 <- qbeta(c(0.005, 0.995), sh, sh)
        return(abs(diff(q99)))
    })
}
# Proportion of population that's reproducing:
p_repro <- function(shape, offset) {
    mapply(\(sh, of) {
        a0 <- beta_starts(sh, of, 1)
        return(sum(a0[10:29]))
    }, sh = shape, of = offset)
}

# What're the "real" values of the Weibull fecundity parameters?
# fop ----
fop <- optim(c(2, 8),
             function(pars, fecunds, max_f, return_fit = FALSE) {
                 if (any(pars < 0)) return(1e10)
                 stages <- 0:length(fecunds)
                 p_distr_vals <- pweibull(stages, shape = pars[1], scale = pars[2])
                 fit_fecunds <- diff(p_distr_vals) / sum(diff(p_distr_vals))
                 fit_fecunds <- fit_fecunds * max_f / max(fit_fecunds)
                 if (return_fit) return(fit_fecunds)
                 return(sum(abs(fit_fecunds - fecunds), na.rm = TRUE))
             },
             fecunds = line_s$leslie[1,9:29,1],
             max_f = max(line_s$leslie[,,1]))


test_sim <- function(i, .known_L, .match_lambda, .fit_survs, .method = "optim") {
    shape <- runif(1, 0.5, 6)
    offset <- runif(1, 0, 1)
    K <- 1800
    sim_df <- sim_aphids(shape, offset, K)
    fits <- fit_sims(sim_df, known_L = .known_L, match_lambda = .match_lambda,
                     fit_survs = .fit_survs)
    if (length(fits) == 1) {
        if (fits == "multiple models") {
            stop(sprintf(paste("Multiple equally well-fitted models with",
                               "different parameter values for",
                               "shape = %s and offset = %s"), shape, offset))
        }
        if (fits == "no converged") {
            stop(sprintf("is.null(op)\n shape = %s\n offset = %s",
                         shape, offset))
        }
        if (fits == "low rows") {
            stop(sprintf("nrow(sim_df) < 3L\n shape = %s\n offset = %s",
                         shape, offset))
        }
        if (fits == "not a model") {
            stop(sprintf("!is.list(op)\n shape = %s\n offset = %s",
                         shape, offset))
        }
        if (fits == "not converged") {
            stop(sprintf("op[[\"convergence\"]] != 0\n shape = %s\n offset = %s",
                         shape, offset))
        }
    }
    if (any(is.na(fits))) {
        stop(sprintf("any(is.na(fits))\n shape = %s\n offset = %s",
                     shape, offset))
    }
    .obs_vars <- c(shape, offset)
    if (!.known_L) .obs_vars <- c(.obs_vars, fop$par)
    if (.fit_survs) .obs_vars <- c(.obs_vars, 0.8329342)
    tibble(obs = .obs_vars,
           fit = unname(fits),
           param = names(fits),
           rep = i) |>
        mutate(param = factor(param, levels = param))
}
test_all_sims <- function(n_sims, .known_L, .match_lambda, .fit_survs) {
    p <- progressor(steps = n_sims)
    future_lapply(1:n_sims, function(i, ...) {
        ifits <- test_sim(i, .known_L, .match_lambda, .fit_survs)
        p()
        return(ifits)
    },
    future.seed = TRUE,
    future.packages = c("tidyverse", "gameofclones", "aphidsync"))
}




# TESTING OTHER OPTIMIZERS -----
shape <- runif(1, 0.5, 6)
offset <- runif(1, 0, 1)
sim_df <- sim_aphids(shape, offset, 1800)
cycles = TRUE
known_L = TRUE
cpp = TRUE
match_lambda = FALSE
fit_survs = TRUE
fit_L_control = list()
# method = "optim"
method = "BOBYQA"
optim_control = list()



fecund <- max(line_s$leslie[,,1])

re_df <- sim_df |>
    arrange(time) |>
    mutate(re = (lead(N) / N)^(1/(lead(time) - time)))

val <- 1e10
op <- NULL
regr_ptr <- NULL
fit_L_args <- list(L_method = "BOBYQA",
                   L_upper_bound = 1000,
                   L_tol = 1e-4,
                   L_max_iters = 1000,
                   max_shape = 800)
optim_control <- list()
optimize_f <- NULL
if (method == "optim") {
    optimize_f <- optim
} else if (method == "BOBYQA") {
    optimize_f <- nloptr::bobyqa
} else if (method == "COBYLA") {
    optimize_f <- nloptr::cobyla
} else if (method == "NELDERMEAD") {
    optimize_f <- nloptr::neldermead
} else if (method == "SBPLX") {
    optimize_f <- nloptr::sbplx
} else stop("unknown method")
if (cpp) {
    regr_ptr <- make_regr_ptr()
}
if (length(fit_L_control) > 0) {
    stopifnot(all(names(fit_L_control) %in% names(fit_L_args)))
    for (n in names(fit_L_control)) fit_L_args[[n]] <- fit_L_control[[n]]
}
if (method != "optim") {
    optim_control[["stopval"]] <- 0 # stop minimization at this value
    optim_control[["xtol_rel"]] <- 1e-8 # stop on small optimization step
    optim_control[["ftol_rel"]] <- 1e-6 # stop on change times function value
}
if (length(optim_control) > 0) {
    for (n in names(optim_control)) optim_control[[n]] <- optim_control[[n]]
}
n_starts <- 100L
if (known_L) {
    known_L_mat <- line_s$leslie[,,1]
} else {
    known_L_mat <- matrix(0, 0, 0)
}

guess <- runif(5, 0, c(10, 1, 10, 20, 1))

if (cpp) {
    args <- list(par = guess, fn = fit_aphids0, control = optim_control,
                         known_L = known_L_mat,
                         re = re_df$re, time = re_df$time,
                         fecund = fecund, regr_ptr = regr_ptr,
                         match_lambda = match_lambda,
                         fit_survs = fit_survs,
                         max_shape = fit_L_args[["max_shape"]],
                         L_method = fit_L_args[["L_method"]],
                         L_upper_bound = fit_L_args[["L_upper_bound"]],
                         L_tol = fit_L_args[["L_tol"]],
                         L_max_iters = fit_L_args[["L_max_iters"]])
} else {
    args <- list(par = guess, fn = fit_aphids0R,
                 control = optim_control,
                 fecund = fecund,
                 match_lambda = match_lambda,
                 fit_survs = fit_survs,
                 known_L = known_L, re_df = re_df,
                 L_upper_bound = fit_L_args[["L_upper_bound"]],
                 max_shape = fit_L_args[["max_shape"]])
}
if (method != "optim") {
    args[["x0"]] <- args[["par"]]
    args[["par"]] <- NULL
    args[["lower"]] <- c(0, 0, 1e-10, 1e-10, 0)[1:length(guess)]
    args[["upper"]] <- c(fit_L_args[["max_shape"]], 1, 100, 100, 1)[1:length(guess)]
}

# args[["lower"]] <- c(0, 0, 1e-10, 1e-10, 0)[1:length(guess)]
# args[["upper"]] <- c(fit_L_args[["max_shape"]], 1, 100, 100, 1)[1:length(guess)]
args[["lower"]] <- NULL
args[["upper"]] <- NULL
args[["control"]] <- list()
args[["control"]][["stopval"]] <- 0 # stop minimization at this value
args[["control"]][["xtol_rel"]] <- 1e-8 # stop on small optimization step
# args[["control"]][["ftol_rel"]] <- 1e-6 # stop on change times function value
args[["control"]][["ftol_abs"]] <- 0
new_op <- do.call(nloptr::sbplx, args)
new_op


fit_aphids0(pars = new_op[["par"]],
known_L = known_L_mat,
re = re_df$re, time = re_df$time,
fecund = fecund, regr_ptr = regr_ptr,
match_lambda = match_lambda,
fit_survs = fit_survs,
max_shape = fit_L_args[["max_shape"]],
L_method = fit_L_args[["L_method"]],
L_upper_bound = fit_L_args[["L_upper_bound"]],
L_tol = fit_L_args[["L_tol"]],
L_max_iters = fit_L_args[["L_max_iters"]])


args2 <- args
args2[["par"]] <- args2[["x0"]]
args2[["x0"]] <- NULL
args2[["lower"]] <- NULL
args2[["upper"]] <- NULL
args2[["control"]] <- NULL
new_op2 <- do.call(optim, args2)


new_op <- do.call(optimize_f, args)

new_op



# $par
# [1]  1.1188251  0.9094286 30.4121541 33.6112852  0.4278546
#
# $value
# [1] 1.221916e-07
#
# $counts
# function gradient
# 462       NA
#
# $convergence
# [1] 0
#
# $message
# NULL

# if (new_op$convergence == 0 && new_op$value < val) {
#     val <- new_op$value
#     op <- new_op
# }
# if (new_op$convergence == 0 && !is.null(op) && new_op$value == val) {
#     if (!isTRUE(all.equal(op$par, new_op$par))) {
#         # This triggers a specific error in `test_sim`:
#         return("multiple models")
#     }
# }









# WEIRDNESS ----
# sh <- 3
# of <- 0.5
# wsh <- fop$par[1]
# wsc <- fop$par[2]
# spa <- 0.8329342
# obs <- c(sh, of, wsh, wsc, spa)
#
# max_f <- max(line_s$leslie[,,1])
#
# sim_df <- sim_aphids(sh, of, 1800)
# re_df <- sim_df |>
#     arrange(time) |>
#     mutate(re = (lead(N) / N)^(1/(lead(time) - time)))
#
# fits <- fit_sims(sim_df, known_L = FALSE, match_lambda = FALSE)
# cbind(fits, obs)
#
# fit_aphids0R(fits, FALSE, re_df, max_f, FALSE)
# fit_aphids0R(obs, FALSE, re_df, max_f, FALSE)









# simulations ----

overwrite <- TRUE

# Takes 25 sec
if (overwrite || ! file.exists("_scripts/known_fit_df.rds")) {
    t0 <- Sys.time()
    set.seed(1740563097)
    # Note that .match_lambda par doesn't matter here.
    known_fit_df <- test_all_sims(100L, .known_L = TRUE, .match_lambda = TRUE,
                                  .fit_survs = FALSE) |>
        do.call(what = bind_rows)
    cat("Finished #1!\n")
    t1 <- Sys.time()
    write_rds(known_fit_df, "_scripts/known_fit_df.rds")
    print(t1 - t0); rm(t0, t1)
} else known_fit_df <- read_rds("_scripts/known_fit_df.rds")


# Takes ~22 min
if (overwrite || ! file.exists("_scripts/unknown_fits_df.rds")) {
    t0 <- Sys.time()
    set.seed(1188441535)
    unknown_fits_df <- crossing(.ml = c(TRUE, FALSE),
                                .fs = c(TRUE, FALSE)) |>
        pmap(\(.ml, .fs) {
            X <- test_all_sims(100L, .known_L = FALSE, .match_lambda = .ml,
                               .fit_survs = .fs) |>
                do.call(what = bind_rows) |>
                mutate(match_lambda = .ml, fit_survs = .fs)
            cat(sprintf("Finished match_lambda=%s, fit_survs=%s!\n", .ml, .fs))
            return(X)
        }) |>
        do.call(what = bind_rows) |>
        mutate(match_lambda = factor(match_lambda, levels = c(TRUE, FALSE),
                                     labels = c("lambda", "max_f")),
               fit_survs = factor(fit_survs, levels = c(TRUE, FALSE),
                                  labels = c("fit_survs", "survs=1")))
    t1 <- Sys.time()
    # write_rds(unknown_fits_df, "_scripts/unknown_fits_df.rds")
    print(t1 - t0); rm(t0, t1)
} else unknown_fits_df <- read_rds("_scripts/unknown_fits_df.rds")


unknown_fits_df0 <- read_rds("_scripts/unknown_fits_df.rds")



# LEFT OFF HERE ----


bind_rows(unknown_fits_df |> mutate(logit = FALSE),
          unknown_fits_df0 |> mutate(logit = TRUE)) |>
    mutate(logit = factor(logit)) |>
    split(~ logit + match_lambda + fit_survs + rep, drop = TRUE) |>
    map_dfr(\(d) {
        sh <- c(d$obs[d$param == "shape"], d$fit[d$param == "shape"])
        of <- c(d$obs[d$param == "offset"], d$fit[d$param == "offset"])
        d |>
            add_row(obs = c(width99(sh[1], of[1]), p_repro(sh[1], of[1])),
                    fit = c(width99(sh[2], of[2]), p_repro(sh[2], of[2])),
                    param = c("width99", "p_repro"),
                    logit = d$logit[1:2],
                    rep = d$rep[1:2],
                    match_lambda = d$match_lambda[1:2],
                    fit_survs = d$fit_survs[1:2])
    }) |>
    mutate(param = factor(param, levels = c("shape", "offset", "wshape",
                                            "wscale", "spar", "width99",
                                            "p_repro"))) |>
    group_by(logit, param, match_lambda, fit_survs) |>
    summarize(ss = mean(abs(obs - fit)), .groups = "drop") |>
    mutate(id = paste(match_lambda, fit_survs) |> factor()) |>
    ggplot(aes(ss, id, color = id)) +
    geom_vline(xintercept = 0, linewidth = 0.5, color = "gray70") +
    geom_segment(aes(xend  = 0, yend = id), linewidth = 1) +
    geom_point(aes(shape = logit), size = 4) +
    facet_wrap(~ param, scales = "free_x") +
    scale_color_viridis_d(option = "magma", begin = 0.1, end = 0.9, guide = "none") +
    scale_shape_manual(values = c(1, 6)) +
    xlab("mean(| obs - fit |)")  +
    ylab("Method")



unknown_fits_df |>
    split(~ match_lambda + fit_survs + rep) |>
    map_dfr(\(d) {
        sh <- c(d$obs[d$param == "shape"], d$fit[d$param == "shape"])
        of <- c(d$obs[d$param == "offset"], d$fit[d$param == "offset"])
        d |>
            add_row(obs = c(width99(sh[1], of[1]), p_repro(sh[1], of[1])),
                    fit = c(width99(sh[2], of[2]), p_repro(sh[2], of[2])),
                    param = c("width99", "p_repro"),
                    rep = d$rep[1:2],
                    match_lambda = d$match_lambda[1:2],
                    fit_survs = d$fit_survs[1:2])
    }) |>
    mutate(param = factor(param, levels = c("shape", "offset", "wshape",
                                            "wscale", "spar", "width99",
                                            "p_repro"))) |>
    group_by(param, match_lambda, fit_survs) |>
    summarize(ss = mean(abs(obs - fit)), .groups = "drop") |>
    mutate(id = paste(match_lambda, fit_survs) |> factor()) |>
    ggplot(aes(ss, id, color = id)) +
    geom_vline(xintercept = 0, linewidth = 0.5, color = "gray70") +
    geom_segment(aes(xend  = 0, yend = id), linewidth = 1) +
    geom_point(size = 4) +
    facet_wrap(~ param, scales = "free_x") +
    scale_color_viridis_d(option = "magma", begin = 0.1, end = 0.9, guide = "none") +
    xlab("mean(| obs - fit |)")  +
    ylab("Method")



unknown11 <- c(levels(unknown_fits_df$param)[1:2], "width99", "p_repro") |>
    set_names() |>
    map(\(.p) {
        if (.p %in% levels(unknown_fits_df$param)) {
            d <- unknown_fits_df |>
                filter(param == .p)
        } else {
            stopifnot(.p %in% c("width99", "p_repro"))
            f <- eval(parse(text = .p))
            d <- unknown_fits_df |>
                filter(param %in% c("shape", "offset")) |>
                pivot_wider(names_from = "param", values_from = c(obs, fit)) |>
                mutate(obs = f(obs_shape, obs_offset),
                       fit = f(fit_shape, fit_offset))
        }
        d |>
            ggplot(aes(obs, fit)) +
            geom_abline(slope = 1, intercept = 0, linetype = 2, color = "red") +
            geom_point() +
            facet_wrap(~ match_lambda + fit_survs, nrow = 2, scales = "free_y") +
            ggtitle(.p) +
            theme(plot.title = element_text(size = 16, hjust = 0.5))
    })

unknown_hist <- levels(unknown_fits_df$param)[3:5] |>
    set_names() |>
    map(\(.p) {
        .obs <- unknown_fits_df$obs[[which(unknown_fits_df$param == .p)[1]]]
        unknown_fits_df |>
            filter(param == .p) |>
            ggplot(aes(fit)) +
            geom_histogram(aes(x = fit), bins = 25, fill = "dodgerblue3") +
            geom_vline(data = tibble(obs = .obs),
                       aes(xintercept = obs), linetype = 2, color = "firebrick1") +
            facet_wrap(~ match_lambda + fit_survs, nrow = 2, scales = "free_y") +
            ggtitle(.p) +
            scale_y_continuous() +
            theme(plot.title = element_text(size = 16, hjust = 0.5))
    })

levels(unknown_fits_df$param)

unknown11[["shape"]]
unknown11[["offset"]]
unknown11[["width99"]]
unknown11[["p_repro"]]


unknown_hist[["wshape"]]
unknown_hist[["wscale"]]
unknown_hist[["spar"]]








# plots ----

known_fit_df |>
    ggplot(aes(obs, fit)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "red") +
    geom_point() +
    # geom_point() +
    facet_wrap(~ param, nrow = 1, scales = "free") +
    ggtitle("with known Leslie matrix") +
    # coord_equal() +
    scale_y_continuous()

lambda_fit_df |>
    filter(param %in% c("shape", "offset")) |>
    # filter(fit < 15) |>
    ggplot(aes(obs, fit)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "red") +
    geom_point() +
    facet_wrap(~ param, nrow = 1, scales = "free") +
    ggtitle("with unknown Leslie matrix (using lambda)") +
    # coord_equal() +
    scale_y_continuous()

max_fit_df |>
    filter(param %in% c("shape", "offset")) |>
    # filter(fit < 15) |>
    ggplot(aes(obs, fit)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "red") +
    geom_point() +
    facet_wrap(~ param, nrow = 1, scales = "free") +
    ggtitle("with unknown Leslie matrix (using max)") +
    # coord_equal() +
    scale_y_continuous()




lambda_fit_df |>
    filter(!param %in% c("shape", "offset")) |>
    # filter(fit < 15) |>
    ggplot(aes(fit)) +
    geom_histogram(aes(x = fit), bins = 25, fill = "dodgerblue3") +
    geom_vline(data = tibble(obs = fop$par,
                             param = levels(lambda_fit_df$param)[-1:-2]) |>
                   mutate(param = factor(param, levels = param)),
               aes(xintercept = obs), linetype = 2, color = "firebrick1") +
    facet_wrap(~ param, nrow = 1, scales = "free") +
    ggtitle("with unknown Leslie matrix (lambda)") +
    # coord_equal() +
    scale_y_continuous()

max_fit_df |>
    filter(!param %in% c("shape", "offset")) |>
    # filter(fit < 15) |>
    ggplot(aes(fit)) +
    geom_histogram(aes(x = fit), bins = 25, fill = "dodgerblue3") +
    geom_vline(data = tibble(obs = fop$par,
                             param = levels(lambda_fit_df$param)[-1:-2]) |>
                   mutate(param = factor(param, levels = param)),
               aes(xintercept = obs), linetype = 2, color = "firebrick1") +
    facet_wrap(~ param, nrow = 1, scales = "free") +
    ggtitle("with unknown Leslie matrix (max)") +
    # coord_equal() +
    scale_y_continuous()



fits <- with(list(shape = median(lambda_fit_df$fit[lambda_fit_df$param == "wshape"]),
                  scale = median(lambda_fit_df$fit[lambda_fit_df$param == "wscale"]),
                  lambda = max(abs(eigen(line_s$leslie[,,1])[["values"]]))),
             {
                 x <- calc_L(shape, scale, lambda)[1,1]
                 L <- make_L1(shape, scale)
                 L[1,9:29] <- L[1,9:29] * x
                 L
             })
optims <- with(list(shape = fop$par[1],
                    scale = fop$par[2],
                    lambda = max(abs(eigen(line_s$leslie[,,1])[["values"]]))),
               {
                   x <- calc_L(shape, scale, lambda)[1,1]
                   L <- make_L1(shape, scale)
                   L[1,9:29] <- L[1,9:29] * x
                   L
               })
fits_max <- with(list(shape = median(max_fit_df$fit[max_fit_df$param == "wshape"]),
                      scale = median(max_fit_df$fit[max_fit_df$param == "wscale"]),
                      max_f = max(line_s$leslie[1,9:29,1])),
               {
                   L <- make_L1(shape, scale)
                   L[1,9:29] <- L[1,9:29] * max_f / max(L[1,9:29])
                   L
               })
optims_max <- with(list(shape = fop$par[1],
                       scale = fop$par[2],
                       max_f = max(line_s$leslie[1,9:29,1])),
               {
                   L <- make_L1(shape, scale)
                   L[1,9:29] <- L[1,9:29] * max_f / max(L[1,9:29])
                   L
               })

{
    print(max(abs(eigen(fits)[["values"]])))
    print(max(abs(eigen(optims)[["values"]])))
    print(max(abs(eigen(fits_max)[["values"]])))
    print(max(abs(eigen(optims_max)[["values"]])))
    print(max(abs(eigen(line_s$leslie[,,1])[["values"]])))
}


tibble(stage = 9:29,
       fit = fits[1,9:29],
       optim = optims[1,9:29],
       fit_max = fits_max[1,9:29],
       optim_max = optims_max[1,9:29],
       real = line_s$leslie[1,9:29,1]) |>
    pivot_longer(-stage, names_to = "source") |>
    # mutate(source = factor(source, levels = c("fit", "optim", "real"))) |>
    ggplot(aes(stage, value, color = source)) +
    geom_point() +
    geom_line() +
    scale_color_viridis_d(begin = 0.2)

plot(line_s$leslie[,,1][row(line_s$leslie[,,1]) - col(line_s$leslie[,,1]) == 1], type = "l", ylab = "survival", xlab = "stage")







for (j in 1:2) {
    f <- list(width99, p_repro)[[j]]
    st <- c("99th quartile width", "Proportion reproducing")[[j]]
    p <- imap(list("known Leslie" = known_fit_df,
              # "unknown Leslie (lambda)" = lambda_fit_df,
              "unknown Leslie (max)" = max_fit_df),
         \(x, i) {
             x |>
                 filter(param %in% c("shape", "offset")) |>
                 pivot_wider(names_from = "param", values_from = c(obs, fit)) |>
                 mutate(obs = f(obs_shape, obs_offset),
                        fit = f(fit_shape, fit_offset)) |>
                 ggplot(aes(obs, fit)) +
                 geom_abline(slope = 1, intercept = 0, linetype = 2, color = "red") +
                 geom_point() +
                 labs(title = i) +
                 # coord_equal() +
                 scale_y_continuous() +
                 theme(plot.title = element_text(size = 12, hjust = 0.5))
         }) |>
        c(list(nrow = 1)) |>
        do.call(what = wrap_plots) +
        plot_annotation(title = st,
                        theme = theme(plot.title = element_text(hjust = 0.5,
                                                                size = 16)))
    plot(p)
}; rm(j, f, st, p)





# survivals ----

surv_df <- tibble(surv = line_s$leslie[,,1][row(line_s$leslie[,,1]) -
                                                col(line_s$leslie[,,1]) == 1],
                  stage = 1:28,
                  stage0 = stage - 28L,
                  surv1 = surv - 1)

# s_mod <- lm(surv1 ~ poly(stage, 2, raw = TRUE) + 0, surv_df)
# s_mod <- lm(surv ~ lhs, surv_df)
s_mod <- lm(surv1 ~ 0 + I(1 / sqrt(1 + stage0^2)), surv_df)


surv_df |>
    ggplot(aes(stage, surv)) +
    geom_hline(yintercept = c(0, 1), color = "gray70") +
    geom_line() +
    geom_line(data = surv_df |>
                  mutate(surv = predict(s_mod) + 1),
                  # mutate(surv = 1 + 0.0003926 * stage0 +
                  #            -0.8278490 /sqrt(1 + stage0^2)),
              color = "dodgerblue", linetype = 2)
