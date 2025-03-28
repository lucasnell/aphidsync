
suppressPackageStartupMessages({
    library(tidyverse)
    library(aphidsync)
    library(future.apply)
    library(progressr)
    library(patchwork)
    library(ggtext)
})

plan(multisession, workers = max(parallel::detectCores() - 2L, 1L))
handlers(global = TRUE)
handlers("progress")


if (interactive() && file.exists(".Rprofile")) source(".Rprofile")


# File created here to store simulation data:
rds_files <- list(with_K = "_data/simple_vary_fs_fit_df.rds",
                  no_K = "_data/simple_vary_fs_noK_fit_df.rds",
                  with_K_adj1 = "_data/simple_vary_fs_adj1_fit_df.rds",
                  no_K_adj1 = "_data/simple_vary_fs_noK_adj1_fit_df.rds",
                  with_K_adj2 = "_data/simple_vary_fs_adj2_fit_df.rds",
                  no_K_adj2 = "_data/simple_vary_fs_noK_adj2_fit_df.rds")


col_pal <- c(resistant = "#218F8DFF",
             susceptible = "#DDE318FF",
             wasp = "#440154FF")

param_lvls <- c("shape", "offset", "K", "width99", "med_age")

# Shared options for `winnowing_optim`:
optim_args__ <- list(box_optim = nloptr::bobyqa,
                     fine_optim = nloptr::bobyqa,
                     polished_optim = nloptr::bobyqa,
                     box_control = list(maxit = 500L, reltol = 1e-6),
                     fine_control = list(maxit = 1000L, reltol = 1e-8),
                     polished_control = list(maxit = 10e3L, reltol = 1e-10),
                     n_fine = 500L,
                     n_polished = 250L)



# Adjust survivals stochastically based on input standard deviation.
# Random variables are (optionally) normalized so that they always have a
# mean of zero.
adjust_survs <- function(L, sigma_s, normalize = FALSE) {
    survs <- L[row(L) - col(L) == 1]
    rnds <- rnorm(length(survs), 0, sigma_s)
    if (normalize) rnds <- rnds - mean(rnds)
    survs <- inv_logit(logit(survs) + rnds)
    L[row(L) - col(L) == 1] <- survs
    return(L)
}
# Same thing for fecundities.
# Here, random variables are (optionally) normalized so that they always have a
# mean of 1 because we're multiplying them by the fecundities, not adding them.
adjust_fecunds <- function(L, sigma_f, normalize = FALSE) {
    rnds <- exp(rnorm(ncol(L), 0, sigma_f))
    if (normalize) rnds <- rnds / mean(rnds)
    L[1,] <- L[1,] * rnds
    return(L)
}


L0 <- matrix(0, 5, 5)
L0[row(L0) - col(L0) == 1] <- 0.96 * exp(0.2 * 0:-3)
L0[1,5] <- 10

L1 <- L0
# Use optimize to get value that closely matches overall growth rate for `L0`:
L1[1,4:5] <- optimize(\(x) {
    L <- L0
    L[1,4:5] <- x
    abs(calc_lambda(L0) - calc_lambda(L))
}, c(1, L0[1,5]))[["minimum"]]

calc_lambda(L0); calc_lambda(L1)






# =============================================================================*
# functions ----
# =============================================================================*

sim_aphids <- function(shape, offset, K, L,
                       sample_filter = FALSE,
                       max_t = 150) {

    aphids0 <- beta_starts(shape, offset, 100, nrow(L))
    .time <- 0:max_t

    sim_df <- tibble(N = as.numeric(sim_N(aphids0, L, .time, K)),
                     time = .time)

    if (sample_filter){
        sim_df <- sim_df |>
            # Filter for sampling twice per week.
            # Starting with the first day, sample that day of the week and
            # one other three days later:
            filter(time %% 7L == 0 | time %% 7L == 3L)
    }

    return(sim_df)

}




# Test 1 simulation fit for a single fecundity and survival combination
test_sim <- function(i, sim_L, fit_L,
                     no_K = FALSE,
                     adjust_L = 0L,
                     fit_max_t = 30L,
                     optim_args = list()) {
    shape <- runif(1, 0.5, 6)
    offset <- runif(1, 0, 1)
    if (no_K) {
        K <- Inf
    } else K <- runif(1, 1000, 5000)
    re_df <- sim_aphids(shape, offset, K, sim_L)
    fits <- fit_known_leslie(re_df, fit_L, N0 = 100, K = K,
                             adjust_L = adjust_L, optim_args = optim_args)
    if (any(is.na(fits))) {
        stop(sprintf("any(is.na(fits))\n shape = %s\n offset = %s",
                     shape, offset))
    }
    fits <- fits[param_lvls]
    tibble(obs = c(shape, offset, K, width99(shape), med_age(shape, offset, ncol(L))),
           fit = unname(fits),
           param = names(fits),
           rep = i) |>
        mutate(param = factor(param, levels = param_lvls))
}


# Test all 100 simulation fits for a single fecundity and survival combination
test_all_sims <- function(n_sims, L,
                          sigma_s = 0,
                          sigma_f = 0,
                          no_K = FALSE,
                          adjust_L = 0L,
                          fit_max_t = 30L,
                          optim_args = list()) {
    p <- progressor(steps = n_sims)
    out <- future_lapply(1:n_sims, function(i) {
        fit_L <- L
        if (sigma_s > 0) fit_L <- adjust_survs(fit_L, sigma_s)
        if (sigma_f > 0) fit_L <- adjust_fecunds(fit_L, sigma_f)
        ifits <- test_sim(i = i, sim_L = L, fit_L = fit_L, no_K = no_K,
                          adjust_L = adjust_L, fit_max_t = fit_max_t,
                          optim_args = optim_args)
        p()
        return(ifits)
    },
    future.seed = TRUE,
    future.packages = c("tidyverse", "aphidsync")) |>
        do.call(what = bind_rows) |>
        mutate(sigma_s = sigma_s, sigma_f = sigma_f)
    return(out)
}



compare_f <- function(nsims, sigma_f) {
    L <- matrix(0, 5, 5)
    out <- map_dfr(1:nsims, \(i) {
        L[1,] <- runif(5, 0, 10)
        L2 <- adjust_fecunds(L, sigma_f)
        return(tibble(rep = i, sigma_f = sigma_f,
                      before = L[1,], after = L2[1,]))
    })
    return(out)
}
compare_s <- function(nsims, sigma_s) {
    L <- matrix(0, 5, 5)
    s_inds <- which(row(L) - col(L) == 1)
    out <- map_dfr(1:nsims, \(i) {
        L[s_inds] <- runif(length(s_inds))
        L2 <- adjust_survs(L, sigma_s)
        return(tibble(rep = i, sigma_s = sigma_s,
                      before = L[s_inds], after = L2[s_inds]))
    })
    return(out)
}


# These two combinations look good for testing:
sigma_s_vals <- c(0.1, 0.2, 0.5)
sigma_f_vals <- c(0.05, 0.1, 0.2)

# sigma_s_vals |>
#     map_dfr(\(.ss) compare_s(100, .ss)) |>
#     mutate(sigma_s = factor(sigma_s)) |>
#     ggplot(aes(before, after)) +
#     geom_point() +
#     geom_abline(slope = 1, intercept = 0, linetype = "22", color ="red") +
#     facet_wrap(~ sigma_s, nrow = 1) +
#     coord_equal()
# sigma_f_vals |>
#     map_dfr(\(.sf) compare_f(100, .sf)) |>
#     mutate(sigma_f = factor(sigma_f)) |>
#     ggplot(aes(before, after)) +
#     geom_point() +
#     geom_abline(slope = 1, intercept = 0, linetype = "22", color ="red") +
#     facet_wrap(~ sigma_f, nrow = 1) +
#     coord_equal()




# =============================================================================*
# simulations ----
# =============================================================================*



overwrite <- TRUE

# This function times simulation and saves file, returning output.
# It only has arguments for things that differ among the 6 combinations below.
run_write_sim <- function(rds_file, seed, no_K, adjust_L) {
    t0 <- Sys.time()
    set.seed(seed)
    out_df <- crossing(.sigma_s = sigma_s_vals,
                       .sigma_f = sigma_f_vals,
                       .repro = c("terminal", "non-terminal")) |>
        mutate(.i = sprintf("%i of %i", 1:n(), n())) |>
        pmap_dfr(\(.sigma_s, .sigma_f, .repro, .i) {
            if (.repro == "terminal") {
                L <- L0
            } else L <- L1
            out <- test_all_sims(60L, L, sigma_s = .sigma_s,
                                 sigma_f = .sigma_f,
                                 no_K = no_K,
                                 adjust_L = adjust_L,
                                 optim_args = optim_args__) |>
                mutate(repro = .repro)
            cat(sprintf(paste("finished sigma_s = %.1f, sigma_f = %.1f,",
                              "repro = %s (%s)\n"),
                        .sigma_s, .sigma_f, .repro, .i))
            return(out)
        })
    write_rds(out_df, rds_file)
    t1 <- Sys.time()
    cat("==================\n")
    cat("TOTAL TIME\n")
    cat(sprintf("%.2f minutes\n", as.numeric(difftime(t1, t0, units = "min"))))
    cat("==================\n")
    return(out_df)
}


# Each of the simulations below take ~9 minutes.


# -----------------*
# No Leslie adjustment
# -----------------*
# With carrying capacity:
if (overwrite || ! file.exists(rds_files$with_K)) {
    fit_test_df <- run_write_sim(rds_files$with_K, 666175554,
                                 no_K = FALSE, adjust_L = 0L)
} else fit_test_df <- read_rds(rds_files$with_K)

# Without carrying capacity:
if (overwrite || ! file.exists(rds_files$no_K)) {
    fit_nok_test_df <- run_write_sim(rds_files$no_K, 1013619437,
                                     no_K = TRUE, adjust_L = 0L)
} else fit_nok_test_df <- read_rds(rds_files$no_K)


# -----------------*
# 1-parameter Leslie adjustment
# -----------------*
# With carrying capacity:
if (overwrite || ! file.exists(rds_files$with_K_adj1)) {
    fit_adj1_test_df <- run_write_sim(rds_files$with_K_adj1, 980091274,
                                      no_K = FALSE, adjust_L = 1L)
} else fit_adj1_test_df <- read_rds(rds_files$with_K_adj1)

# Without carrying capacity:
if (overwrite || ! file.exists(rds_files$no_K_adj1)) {
    fit_nok_adj1_test_df <- run_write_sim(rds_files$no_K_adj1, 832073052,
                                     no_K = TRUE, adjust_L = 1L)
} else fit_nok_adj1_test_df <- read_rds(rds_files$no_K_adj1)


# -----------------*
# 2-parameter Leslie adjustment
# -----------------*
# With carrying capacity:
if (overwrite || ! file.exists(rds_files$with_K_adj2)) {
    fit_adj2_test_df <- run_write_sim(rds_files$with_K_adj2, 236120720,
                                      no_K = FALSE, adjust_L = 2L)
} else fit_adj2_test_df <- read_rds(rds_files$with_K_adj2)

# Without carrying capacity:
if (overwrite || ! file.exists(rds_files$no_K_adj2)) {
    fit_nok_adj2_test_df <- run_write_sim(rds_files$no_K_adj2, 614944245,
                                          no_K = TRUE, adjust_L = 2L)
} else fit_nok_adj2_test_df <- read_rds(rds_files$no_K_adj2)






# =============================================================================*
# plots ----
# =============================================================================*

#' Plotter for varying fecundities and survivals.
#'
#' @param .p Single string indicating parameter to plot.
#' @param .r Single string indicating the reproduction type (`"terminal"` or
#'     `"non-terminal"`) to plot.
#' @param .df Dataframe containing data to plot.
#' @param .ts Suffix to add to title. Defaults to `""`.
#'
vary_fs_plotter <- function(.p, .r, .df, .ts = "") {

    stopifnot(.p %in% param_lvls)
    stopifnot(.p %in% .df$param)
    stopifnot(.r %in% .df$repro)

    g <- function(x, n) factor(x, levels = sort(unique(x)),
                               labels = paste(c("low", "mid", "high"), n))
    pm <- list(shape = "Initial Beta shape",
               offset = "Initial Beta offset",
               K = "Density dependence",
               width99 = "Width of 99% quantile",
               med_age = "Median starting age")

    d <- .df |>
        filter(repro == .r, param == .p)
    .m <- min(d$obs, d$fit)
    p <- d |>
        mutate(sigma_s = g(sigma_s, "<i>&sigma;<sub>s</sub></i>"),
               sigma_f = g(sigma_f, "<i>&sigma;<sub>f</sub></i>")) |>
        ggplot(aes(obs, fit)) +
        # geom_hline(yintercept = .m, color = "gray70", linewidth = 0.5) +
        # geom_vline(xintercept = .m, color = "gray70", linewidth = 0.5) +
        geom_abline(slope = 1, intercept = 0, linetype = 2, color = "red") +
        geom_point() +
        xlab("Actual value") +
        ylab("Fitted value") +
        facet_wrap(~ sigma_s + sigma_f, nrow = 3, scales = "free") +
        ggtitle(sprintf("%s (%s repro)%s", pm[[.p]], .r, .ts)) +
        theme(strip.text = element_markdown())
    return(p)
}





prm <- "width99"

# ----------*
# No Leslie scaling:
# terminal repro:
vary_fs_plotter(prm, "terminal", fit_test_df, ", with K")
vary_fs_plotter(prm, "terminal", fit_nok_test_df, ", no K")
# non-terminal repro:
vary_fs_plotter(prm, "non-terminal", fit_test_df, ", with K")
vary_fs_plotter(prm, "non-terminal", fit_nok_test_df, ", no K")

# ----------*
# 1-param Leslie scaling:
# terminal repro:
vary_fs_plotter(prm, "terminal", fit_adj1_test_df, ", with K (1-param scaling)")
vary_fs_plotter(prm, "terminal", fit_nok_adj1_test_df, ", no K (1-param scaling)")
# non-terminal repro:
vary_fs_plotter(prm, "non-terminal", fit_adj1_test_df, ", with K (1-param scaling)")
vary_fs_plotter(prm, "non-terminal", fit_nok_adj1_test_df, ", no K (1-param scaling)")

# ----------*
# 2-param Leslie scaling:
# terminal repro:
vary_fs_plotter(prm, "terminal", fit_adj2_test_df, ", with K (2-param scaling)")
vary_fs_plotter(prm, "terminal", fit_nok_adj2_test_df, ", no K (2-param scaling)")
# non-terminal repro:
vary_fs_plotter(prm, "non-terminal", fit_adj2_test_df, ", with K (2-param scaling)")
vary_fs_plotter(prm, "non-terminal", fit_nok_adj2_test_df, ", no K (2-param scaling)")











cor_df <- bind_rows(fit_test_df |> mutate(with_k = "with K", scaling = 0L),
          fit_nok_test_df |> mutate(with_k = "no K", scaling = 0L),
          fit_adj1_test_df |> mutate(with_k = "with K", scaling = 1L),
          fit_nok_adj1_test_df |> mutate(with_k = "no K", scaling = 1L),
          fit_adj2_test_df |> mutate(with_k = "with K", scaling = 2L),
          fit_nok_adj2_test_df |> mutate(with_k = "no K", scaling = 2L)) |>
    mutate(scaling = factor(scaling, levels = 0:2,
                            labels = c("no scaling",
                                       paste0(1:2, "-param scaling")))) |>
    # mutate(sigma_s = factor(sigma_s, levels = sort(unique(sigma_s)),
    #                         labels = paste(c("low", "mid", "high"),
    #                                        "<i>&sigma;<sub>s</sub></i>")),
    #        sigma_f = factor(sigma_f, levels = sort(unique(sigma_f)),
    #                         labels = paste(c("low", "mid", "high"),
    #                                        "<i>&sigma;<sub>f</sub></i>"))) |>
    filter(param != "K") |>
    group_by(with_k, scaling, repro, sigma_s, sigma_f, param) |>
    summarize(r = cor(obs, fit), .groups = "drop")

prm = "med_age"
sgm = "sigma_s"

other_sgm <- ifelse(sgm == "sigma_s", "sigma_f", "sigma_s")

cor_df |>
    filter(param == prm) |>
    mutate(!!sgm := factor(.data[[sgm]]),
           !!other_sgm := factor(.data[[other_sgm]],
                                 levels = sort(unique(.data[[other_sgm]])),
                                 labels = sprintf("<i>&sigma;<sub>%s</sub></i> = %.2f",
                                                  str_remove(other_sgm, "sigma_"),
                                                  sort(unique(.data[[other_sgm]]))))) |>
    ggplot(aes(.data[[sgm]], r)) +
    geom_point(aes(color = scaling),
               position = position_dodge(0.5)) +
    geom_linerange(aes(color = scaling, ymin = 0, ymax = r),
                 position = position_dodge2(0.5)) +
    facet_grid(.data[[other_sgm]] ~ repro + with_k) +
    scale_color_viridis_d(begin = 0.2, end = 0.9) +
    xlab(sprintf("<i>&sigma;<sub>%s</sub></i>", str_remove(sgm, "sigma_"))) +
    scale_y_continuous("*r*", limits = c(0,1)) +
    scale_shape_manual(values = c(16, 17)) +
    scale_linetype_manual(values = c("solid", "22")) +
    theme(axis.title = element_markdown(),
          strip.text.y = element_markdown(angle = 0))




bind_rows(fit_test_df |> mutate(with_k = "with K"),
          fit_nok_test_df |> mutate(with_k = "no K")) |>
    mutate(surv_x = factor(surv_x, levels = sort(unique(surv_x)),
                           labels = paste(c("low", "perfect", "high"), "survival")),
           fecund_x = factor(fecund_x, levels = sort(unique(fecund_x)),
                             labels = paste(c("low", "perfect", "high"), "fecundity"))) |>
    pivot_wider(names_from = "param", values_from = c(obs, fit)) |>
    mutate(obs = med_age(obs_shape, obs_offset),
           fit = med_age(fit_shape, fit_offset)) |>
    select(rep:with_k, obs, fit) |>
    group_by(with_k, repro, surv_x, fecund_x) |>
    summarize(r = cor(obs, fit), .groups = "drop")

