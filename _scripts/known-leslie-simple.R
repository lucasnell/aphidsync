source("_scripts/known-leslie-simple-preamble.R")


# File created here to store simulation data:
rds_files <- list(with_K = "_data/simple_vary_fs_fit_df.rds",
                  no_K = "_data/simple_vary_fs_noK_fit_df.rds",
                  with_K_adj1 = "_data/simple_vary_fs_adj1_fit_df.rds",
                  no_K_adj1 = "_data/simple_vary_fs_noK_adj1_fit_df.rds",
                  with_K_adj2 = "_data/simple_vary_fs_adj2_fit_df.rds",
                  no_K_adj2 = "_data/simple_vary_fs_noK_adj2_fit_df.rds")



# # NOT SURE THESE ARE NECESSARY:
# col_pal <- c(resistant = "#218F8DFF",
#              susceptible = "#DDE318FF",
#              wasp = "#440154FF")
#
# # Shared options for `winnowing_optim`:
# optim_args__ <- list(box_optim = nloptr::bobyqa,
#                      fine_optim = nloptr::bobyqa,
#                      polished_optim = stats::optim,
#                      box_control = list(maxit = 500L, reltol = 1e-6),
#                      fine_control = list(maxit = 1000L, reltol = 1e-8),
#                      polished_control = list(maxit = 10e3L, reltol = 1e-10),
#                      n_fine = 500L,
#                      n_polished = 250L)



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
                                 optim_args = optim_args__,
                                 fit_args = list()) |>
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


# LEFT OFF ----

# TO DO ----
# WRITE SIMULATIONS TO SEE WHAT OPTIMIZATION PARAMETERS RESULT IN THE BEST FITS


#' Seeing what parameters in `winnowing_optim` help avoiding fitting issues.
#'
#' Thus far, I've found that the following doesn't help:
#'   - having stricter control for early fits (box and fine stages)
#'   - having more polished fits
#'   - in isolation, using `bobyqa` instead of `optim` (even when including
#'     lower and upper bounds)
#'   -
#'
#' The following does seem to help:
#'   - having more initial box evaluations
#'   - having more boxes, using `bobyqa` with bounds to speed this up
#'

# plot histograms of differences among top fits:
plot_fits <- function(x, .title = waiver()) {
    x |>
        filter(param != "K") |>
        split(~ param, drop = TRUE) |>
        map(\(x) {
            x |>
                ggplot(aes(fit_diffs)) +
                geom_histogram(aes(y = after_stat(density)), bins = 30) +
                xlab(x$param[[1]])
        }) |>
        c(list(nrow = 2)) |>
        do.call(what = wrap_plots) +
        plot_annotation(title = .title)
}

adjust_L <- 1L
no_K <- FALSE
i <- 4L
n_sims <- 300L

VERY_SMALL <- max(.Machine$double.eps, .Machine$double.neg.eps)
lower_bounds <- c(rep(VERY_SMALL, 3) + c(1,0,0), -formals(fit_known_leslie)[["max_surv_x"]]) |>
    base::`[`(1:(2L+adjust_L)) |>
    aphidsync:::trans_pars(trans_base = formals(fit_known_leslie)[["trans_base"]])
upper_bounds <- c(100, 1-VERY_SMALL,
                  formals(fit_known_leslie)[[ifelse(adjust_L == 1L,
                                                    "max_leslie_x", "max_fecund_x")]],
                  formals(fit_known_leslie)[["max_surv_x"]]) |>
    base::`[`(1:(2L+adjust_L)) |>
    aphidsync:::trans_pars(trans_base = formals(fit_known_leslie)[["trans_base"]])

# With all (or mostly) defaults:
set.seed(13456789); out0 <- test_all_sims(n_sims, LL[[i]],
                                          sigma_s = sigma_s_vals[[2]],
                                          sigma_f = sigma_f_vals[[2]],
                                          no_K = no_K,
                                          adjust_L = adjust_L,
                                          fit_args = list(max_shape = 100))

{
    # Takes ~1 min
    t0 <- Sys.time()
    out <- test_all_sims(n_sims, LL[[i]], sigma_s = sigma_s_vals[[2]],
                         sigma_f = sigma_f_vals[[2]],
                         no_K = no_K,
                         adjust_L = adjust_L,
                         optim_args = list(box_optim = nloptr::bobyqa,
                                           # fine_optim = nloptr::bobyqa,
                                           # polished_optim = nloptr::bobyqa,
                                           # box_control = list(maxit = 500L, reltol = 1e-6),
                                           # fine_control = list(maxit = 1000L, reltol = 1e-8),
                                           # polished_control = list(maxit = 10e3L, reltol = 1e-10),
                                           box_control = list(maxit = 100L, reltol = 1e-4,
                                                              lower = lower_bounds,
                                                              upper = upper_bounds),
                                           # n_bevals = 100L * 10L,
                                           n_boxes = 1000L * 10L,
                                           # n_fine = 100L * 5L,
                                           # n_polished = 20L * 2L,
                                           n_outputs = 3L), # << default value
                         fit_args = list(max_shape = 100))
    t1 <- Sys.time()
    print(t1 - t0); rm(t0, t1)
}


# out0 |> plot_fits("defaults")
# out |> plot_fits("extra stuff")

full_join(out0 |>
              group_by(param) |>
              summarize(out0 = mean(abs(fit_diffs))),
          out |>
              group_by(param) |>
              summarize(out = mean(abs(fit_diffs))),
          by = "param")
 # A tibble: 5 × 3
#   param      out0     out
#   <fct>     <dbl>   <dbl>
# 1 shape   0.132   0.0675
# 2 offset  0.00306 0.00271
# 3 K       0       0
# 4 width99 0.00455 0.00313
# 5 med_age 0.0100  0.00333

full_join(out0 |>
              group_by(param) |>
              summarize(out0 = mean(abs(obs - fit))),
          out |>
              group_by(param) |>
              summarize(out = mean(abs(obs - fit))),
          by = "param")
# # A tibble: 5 × 3
#   param     out0    out
#   <fct>    <dbl>  <dbl>
# 1 shape   0.809  0.775
# 2 offset  0.0374 0.0490
# 3 K       0      0
# 4 width99 0.0225 0.0233
# 5 med_age 0.0667 0.0433





# Using optim for box_optim not as good as bobyqa.
# Using bobyqa (with bounds) and 10e3L boxes (took 13.57967 min):
# # A tibble: 5 × 3
#   param      out0        out
#   <fct>     <dbl>      <dbl>
# 1 shape   0.0622  0.00210
# 2 offset  0.00286 0.00000653
# 3 K       0       0
# 4 width99 0.00335 0.0000263
# 5 med_age 0.00667 0


p1 <- out |>
    filter(param != "med_age") |>
    filter(is.finite(obs)) |>
    ggplot(aes(obs, fit)) +
    # geom_hline(yintercept = .m, color = "gray70", linewidth = 0.5) +
    # geom_vline(xintercept = .m, color = "gray70", linewidth = 0.5) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "red") +
    xlab("Actual value") +
    ylab("Fitted value") +
    facet_wrap(~ param, nrow = 2, scales = "free")
p2 <- out |>
    filter(param == "med_age") |>
    ggplot(aes((obs - fit))) +
    geom_histogram(binwidth = 1, fill = "dodgerblue") +
    geom_vline(xintercept = 0, linetype = 2, color = "red") +
    xlab("Fitted - actual med_age")
p1 + p2 +
    plot_annotation(title = sprintf("%i-stage repro, %i-param scaling, %s", i,
                                    adjust_L, ifelse(no_K, "no K", "with K")))




# Each of the simulations below take ~9 minutes.


# -----------------*
# No Leslie adjustment
# -----------------*
# With carrying capacity:
if (overwrite || ! file.exists(rds_files$with_K)) {
    fit_test_df <- run_write_sim(rds_files$with_K, 666175554,
                                 no_K = FALSE, adjust_L = 0L)
} else fit_test_df <- read_rds(rds_files$with_K)

# # Without carrying capacity:
# if (overwrite || ! file.exists(rds_files$no_K)) {
#     fit_nok_test_df <- run_write_sim(rds_files$no_K, 1013619437,
#                                      no_K = TRUE, adjust_L = 0L)
# } else fit_nok_test_df <- read_rds(rds_files$no_K)


# -----------------*
# 1-parameter Leslie adjustment
# -----------------*
# With carrying capacity:
if (overwrite || ! file.exists(rds_files$with_K_adj1)) {
    fit_adj1_test_df <- run_write_sim(rds_files$with_K_adj1, 980091274,
                                      no_K = FALSE, adjust_L = 1L)
} else fit_adj1_test_df <- read_rds(rds_files$with_K_adj1)

# # Without carrying capacity:
# if (overwrite || ! file.exists(rds_files$no_K_adj1)) {
#     fit_nok_adj1_test_df <- run_write_sim(rds_files$no_K_adj1, 832073052,
#                                      no_K = TRUE, adjust_L = 1L)
# } else fit_nok_adj1_test_df <- read_rds(rds_files$no_K_adj1)


# -----------------*
# 2-parameter Leslie adjustment
# -----------------*
# With carrying capacity:
if (overwrite || ! file.exists(rds_files$with_K_adj2)) {
    fit_adj2_test_df <- run_write_sim(rds_files$with_K_adj2, 236120720,
                                      no_K = FALSE, adjust_L = 2L)
} else fit_adj2_test_df <- read_rds(rds_files$with_K_adj2)

# # Without carrying capacity:
# if (overwrite || ! file.exists(rds_files$no_K_adj2)) {
#     fit_nok_adj2_test_df <- run_write_sim(rds_files$no_K_adj2, 614944245,
#                                           no_K = TRUE, adjust_L = 2L)
# } else fit_nok_adj2_test_df <- read_rds(rds_files$no_K_adj2)






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





prm <- "shape"

# ----------*
# No Leslie scaling:
# terminal repro:
vary_fs_plotter(prm, "terminal", fit_test_df, ", with K")
### vary_fs_plotter(prm, "terminal", fit_nok_test_df, ", no K")
# non-terminal repro:
vary_fs_plotter(prm, "non-terminal", fit_test_df, ", with K")
### vary_fs_plotter(prm, "non-terminal", fit_nok_test_df, ", no K")

# ----------*
# 1-param Leslie scaling:
# terminal repro:
vary_fs_plotter(prm, "terminal", fit_adj1_test_df, ", with K (1-param scaling)")
### vary_fs_plotter(prm, "terminal", fit_nok_adj1_test_df, ", no K (1-param scaling)")
# non-terminal repro:
vary_fs_plotter(prm, "non-terminal", fit_adj1_test_df, ", with K (1-param scaling)")
### vary_fs_plotter(prm, "non-terminal", fit_nok_adj1_test_df, ", no K (1-param scaling)")

# ----------*
# 2-param Leslie scaling:
# terminal repro:
vary_fs_plotter(prm, "terminal", fit_adj2_test_df, ", with K (2-param scaling)")
### vary_fs_plotter(prm, "terminal", fit_nok_adj2_test_df, ", no K (2-param scaling)")
# non-terminal repro:
vary_fs_plotter(prm, "non-terminal", fit_adj2_test_df, ", with K (2-param scaling)")
### vary_fs_plotter(prm, "non-terminal", fit_nok_adj2_test_df, ", no K (2-param scaling)")











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

