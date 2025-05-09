
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


if (interactive() && file.exists(".Rprofile")) source(".Rprofile")

# File created here to store simulation data:
rds_file <- "_data/unknown_fits_df.rds"



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
# Fit the two shape parameters for a set of simulations.
# NOTE: This version doesn't allow for multiple lines (yet).
#
fit_sims <- function(sim_df, fit_survs = TRUE,
                     max_shape = 800) {

    stopifnot(length(unique(sim_df$line)) == 1L)

    if (nrow(sim_df) < 3L) return("low rows") # special meaning in `test_sim`

    # `max_f` parameter is max fecundity
    max_f <- max(line_s$leslie[,,1])

    re_df <- sim_df |>
        arrange(time) |>
        mutate(re = (lead(N) / N)^(1/(lead(time) - time)))

    # Fit across multiple starting parameter values, then choose best-fitting one:
    val <- 1e10
    op <- NULL
    n_starts <- 100L
    if (fit_survs) {
        guess_list <- lapply(1:n_starts, \(i) runif(5, 0, c(10, 1, 10, 20, 1)))
    } else {
        guess_list <- lapply(1:n_starts, \(i) runif(4, 0, c(10, 1, 10, 20)))
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
    if (fit_survs) {
        spar_inds <- sample.int(n_starts, n_hilo)
        for (i in spar_inds[1:(n_hilo %/% 2L)])
            guess_list[[i]][5] <- runif(1, 0, 1e-6)
        for (i in spar_inds[(n_hilo %/% 2L + 1L):n_hilo])
            guess_list[[i]][5] <- runif(1, 1 - 1e-6, 1)
    }
    for (s in 1:n_starts) {
        guess <- guess_list[[s]]
        new_op <- optim(guess, unknown_fit_aphids0,
                        re = re_df$re, time = re_df$time,
                        max_f = max_f,
                        fit_survs = fit_survs,
                        max_shape = max_shape)
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
# $par
# [1] 1.981920 8.254428
#
# $value
# [1] 7.929942
#
# $counts
# function gradient
#       93       NA
#
# $convergence
# [1] 0
#
# $message
# NULL


test_sim <- function(i, .fit_survs) {
    shape <- runif(1, 0.5, 6)
    offset <- runif(1, 0, 1)
    K <- 1800
    sim_df <- sim_aphids(shape, offset, K)
    fits <- fit_sims(sim_df, fit_survs = .fit_survs)
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
    .obs_vars <- c(shape, offset, fop$par)
    if (.fit_survs) .obs_vars <- c(.obs_vars, 0.8329342)
    tibble(obs = .obs_vars,
           fit = unname(fits),
           param = names(fits),
           rep = i) |>
        mutate(param = factor(param, levels = param))
}
test_all_sims <- function(n_sims, .fit_survs) {
    p <- progressor(steps = n_sims)
    future_lapply(1:n_sims, function(i, ...) {
        ifits <- test_sim(i, .fit_survs)
        p()
        return(ifits)
    },
    future.seed = TRUE,
    future.packages = c("tidyverse", "gameofclones", "aphidsync"))
}








overwrite <- FALSE

# simulations ----

# Takes ~4 min
if (overwrite || ! file.exists(rds_file)) {
    t0 <- Sys.time()
    set.seed(1188441535)
    unknown_fits_df <- c(TRUE, FALSE) |>
        map_dfr(\(.fs) {
            X <- test_all_sims(100L, .fit_survs = .fs) |>
                do.call(what = bind_rows) |>
                mutate(fit_survs = .fs)
            cat(sprintf("Finished fit_survs=%s!\n", .fs))
            return(X)
        }) |>
        mutate(fit_survs = factor(fit_survs, levels = c(TRUE, FALSE),
                                  labels = c("fit_survs", "survs=1")))
    t1 <- Sys.time()
    write_rds(unknown_fits_df, rds_file)
    print(t1 - t0); rm(t0, t1)
} else unknown_fits_df <- read_rds(rds_file)







# plots ----



# Map parameter names to more descriptive versions:
param_map <- list(shape = "Initial Beta shape",
                  offset = "Initial Beta offset",
                  wshape = "Fecundity Weibull shape",
                  wscale = "Fecundity Weibull scale",
                  spar = "Survival parameter",
                  width99 = "Width of 99% quantile",
                  p_repro = "Proportion of initial aphids reproducing")
# Use map above on a vector and optionally break lines:
pretty_params <- function(uglies, nc = NULL) {
    pretties <- map_chr(uglies, \(u) param_map[[u]])
    if (!is.null(nc)) {
        stopifnot(is.numeric(nc) && length(nc) == 1 && nc > 0)
        pretties <- map_chr(pretties,
                            \(p) paste(strwrap(p, width = nc), collapse = "\n"))
    }
    return(pretties)
}



unknown_ss_p <- unknown_fits_df |>
    split(~ fit_survs + rep) |>
    map_dfr(\(d) {
        sh <- c(d$obs[d$param == "shape"], d$fit[d$param == "shape"])
        of <- c(d$obs[d$param == "offset"], d$fit[d$param == "offset"])
        d |>
            add_row(obs = c(width99(sh[1], of[1]), p_repro(sh[1], of[1])),
                    fit = c(width99(sh[2], of[2]), p_repro(sh[2], of[2])),
                    param = c("width99", "p_repro"),
                    rep = d$rep[1:2],
                    fit_survs = d$fit_survs[1:2])
    }) |>
    mutate(param = factor(param, levels = c("shape", "offset", "wshape",
                                            "wscale", "spar", "width99",
                                            "p_repro"),
                          labels = pretty_params(c("shape", "offset", "wshape",
                                                   "wscale", "spar", "width99",
                                                   "p_repro"),
                                                 nc = 30L))) |>
    group_by(param, fit_survs) |>
    summarize(ss = mean(abs(obs - fit)), .groups = "drop") |>
    ggplot(aes(ss, fit_survs, color = fit_survs)) +
    geom_vline(xintercept = 0, linewidth = 0.5, color = "gray70") +
    geom_segment(aes(xend  = 0, yend = fit_survs), linewidth = 1) +
    geom_point(size = 4) +
    facet_wrap(~ param, scales = "free_x") +
    scale_color_viridis_d(option = "magma", begin = 0.1, end = 0.9, guide = "none") +
    xlab("mean( | obs - fit | )") +
    ylab("Method")
unknown_ss_p
# ggsave("_figures/unknown_ss.pdf", unknown_ss_p, width = 7, height = 6)


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
            facet_wrap(~ fit_survs, nrow = 1, scales = "free_y") +
            ggtitle(param_map[[.p]]) +
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
            facet_wrap(~ fit_survs, nrow = 2, scales = "free_y") +
            ggtitle(param_map[[.p]]) +
            scale_y_continuous() +
            theme(plot.title = element_text(size = 16, hjust = 0.5))
    })

levels(unknown_fits_df$param)

# for (p in names(unknown11)) {
#     ggsave(sprintf("_figures/1to1-%s.pdf", p), unknown11[[p]], width = 6, height = 6)
# }; rm(p)
# for (p in names(unknown_hist)) {
#     ggsave(sprintf("_figures/hists-%s.pdf", p), unknown_hist[[p]], height = 4,
#            width = ifelse(p == "spar", 3, 6))
# }; rm(p)

unknown11[["shape"]]
unknown11[["offset"]]
unknown11[["width99"]]
unknown11[["p_repro"]]


unknown_hist[["wshape"]]
unknown_hist[["wscale"]]
unknown_hist[["spar"]]



