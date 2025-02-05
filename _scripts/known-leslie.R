
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

# Resistant line: high resistance, low population growth rate
line_r <- clonal_line("resistant",
                      density_0 = 32,
                      resistant = TRUE,
                      surv_paras = 0.57,
                      surv_juv_apterous = "low",
                      surv_adult_apterous = "low",
                      repro_apterous = "low",
                      p_instar_smooth = 0)
for(i in 1:3) line_r$leslie[1,,i] <- 0.8 * line_r$leslie[1,,i]









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



fit_LR <- function(pars, max_f, fit_survs) {
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

    x <- max_f / max(L[1, adult_stage:n_stages])

    L[1, adult_stage:n_stages] <- L[1, adult_stage:n_stages] * x

    return(L)
}






#
# Use within optim to find shape parameters that minimize the differences
# between observed and predicted Re values.
#
fit_aphids0R <- function(pars, known_L_mat, re_df, max_f, fit_survs,
                         max_shape = 800) {

    # pars <- c(533.764, 393.973)

    .K <- 1800

    if (any(pars < 0) || pars[2] > 1) return(1e10)
    .shape <- pars[[1]]
    .offset <- pars[[2]]
    # .K <- pars[[3]]
    if (.shape > max_shape) return(1e10)

    aphids0 <- beta_starts(.shape, .offset, total_aphids0 = 32)

    known_L <- is.matrix(known_L_mat) && identical(dim(known_L_mat), c(29L,29L))

    if (known_L) {
        L <- known_L_mat
    } else {
        L <- fit_LR(pars, max_f, fit_survs)
        # `L` is 1e10 if something went wrong in fit_LR
        if (length(L) == 1) return(L)
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

    sae <- sum(abs(re_pred - re_df$re), na.rm = TRUE)
    if (is.na(sae) | sae > 1e10) sae <- 1e10

    return(sae)

}







#
# Fit the two shape parameters for a set of simulations.
# NOTE: This version doesn't allow for multiple lines (yet).
#
fit_sims <- function(sim_df, known_L_mat,
                     fit_survs = TRUE, cpp = TRUE,
                     max_shape = 800, cycles = TRUE) {

    stopifnot(length(unique(sim_df$line)) == 1L)

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

    # `max_f` parameter is max fecundity
    max_f <- max(line_s$leslie[,,1])


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
    n_starts <- 100L
    known_L <- is.matrix(known_L_mat) && identical(dim(known_L_mat), c(29L,29L))
    if (known_L) {
        guess_list <- lapply(1:n_starts, \(i) runif(2, 0, c(10, 1)))
    } else {
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
            new_op <- optim(guess, fit_aphids0,
                            known_L_mat = known_L_mat,
                            re = re_df$re, time = re_df$time,
                            max_f = max_f,
                            fit_survs = fit_survs,
                            max_shape = max_shape)
        } else {
            new_op <- optim(guess, fit_aphids0R, max_f = max_f,
                            fit_survs = fit_survs,
                            known_L_mat = known_L_mat, re_df = re_df,
                            max_shape = max_shape)
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
              color = "dodgerblue", linetype = 2)


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


test_sim <- function(i, .known_L_mat, .fit_survs) {
    shape <- runif(1, 0.5, 6)
    offset <- runif(1, 0, 1)
    K <- 1800
    sim_df <- sim_aphids(shape, offset, K)
    fits <- fit_sims(sim_df, known_L_mat = .known_L_mat, fit_survs = .fit_survs)
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
    .known_L <- is.matrix(.known_L_mat) && identical(dim(.known_L_mat), c(29L,29L))
    .obs_vars <- c(shape, offset)
    if (!.known_L) .obs_vars <- c(.obs_vars, fop$par)
    if (.fit_survs) .obs_vars <- c(.obs_vars, 0.8329342)
    tibble(obs = .obs_vars,
           fit = unname(fits),
           param = names(fits),
           rep = i) |>
        mutate(param = factor(param, levels = param))
}
test_all_sims <- function(n_sims, .known_L_mat, .fit_survs) {
    p <- progressor(steps = n_sims)
    future_lapply(1:n_sims, function(i, ...) {
        ifits <- test_sim(i, .known_L_mat, .fit_survs)
        p()
        return(ifits)
    },
    future.seed = TRUE,
    future.packages = c("tidyverse", "gameofclones", "aphidsync"))
}




# simulations ----

overwrite <- TRUE

# Takes 52 sec
if (overwrite || ! file.exists("_scripts/known_fit_df.rds")) {
    t0 <- Sys.time()
    set.seed(1740563097)
    # Note that .match_lambda par doesn't matter here.
    known_fit_df <- test_all_sims(100L, .known_L_mat = line_s$leslie[,,1],
                                  .fit_survs = FALSE) |>
        do.call(what = bind_rows)
    cat("Finished #1!\n")
    t1 <- Sys.time()
    write_rds(known_fit_df, "_scripts/known_fit_df.rds")
    print(t1 - t0); rm(t0, t1)
} else known_fit_df <- read_rds("_scripts/known_fit_df.rds")





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


