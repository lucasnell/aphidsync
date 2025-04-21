
suppressPackageStartupMessages({
    library(tidyverse)
    library(gameofclones)
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


# Files created here to store simulation data:
rds_files <- list(known = "_data/known_fit_df.rds",
                  vary_fs = "_data/vary_fs_fit_df.rds",
                  vary_fs2 = "_data/vary_fs_fit2_df.rds",
                  vary_fs_samp34 = "_data/vary_fs_samp34_fit_df.rds",
                  vary_fs_term = "_data/vary_fs_term_fit_df.rds",
                  vary_fs_knownK = "_data/vary_fs_knownK_fit_df.rds",
                  vary_fs_N = "_data/vary_fs_N_fit_df.rds")


col_pal <- c(resistant = "#218F8DFF",
             susceptible = "#DDE318FF",
             wasp = "#440154FF")


param_lvls <- c("shape", "offset", "K", "width99", "med_age")

# These two combinations look good for testing:
sigma_s_vals <- c(0.1, 0.2, 0.5)
sigma_f_vals <- c(0.05, 0.1, 0.2)

# Shared options for `winnowing_optim`:
optim_args__ <- list(n_bevals = 1000L)

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







sim_aphids <- function(shape, offset, K, line,
                       sample_filter = FALSE,
                       sd_error = 0,
                       .max_t = 250) {

    line$density_0[,1] <- beta_starts(shape, offset, sum(line$density_0),
                                       nrow(line$density_0))
    line$density_0[,2] <- 0.0

    sim_df <- gameofclones:::sim_gameofclones_full(line,
                                                   n_fields = 1,
                                                   # n_plants = 16,
                                                   n_plants = 1,
                                                   max_t = .max_t,
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
                                                   mean_K = K,  # 1800,
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

    if (sd_error > 0) {
        sim_df$N <- sim_df$N * exp(rnorm(nrow(sim_df), 0, sd_error))
    }

    return(sim_df)

}



# Test 1 simulation fit for a single fecundity and survival combination
test_sim <- function(i, line, fit_L,
                     no_K = FALSE,
                     known_K = TRUE,
                     adjust_L = 1L,
                     fit_max_t = 30L,
                     optim_args = list(),
                     fit_args = list(),
                     sim_args = list()) {
    shape <- runif(1, 1, 6)
    offset <- runif(1, 0, 1)
    if (no_K) {
        K <- Inf
    } else K <- runif(1, 1000, 5000)
    sp_args <- list("shape" = shape, "offset" = offset, "K" = K, "line" = line)
    if (length(sim_args) > 0) for (n in names(sim_args)) sp_args[[n]] <- sim_args[[n]]
    sim_df <- do.call(sim_aphids, sp_args)
    fit_K <- NA_real_
    if (known_K) fit_K <- K
    fkl_args <- list(sim_df = sim_df, L = fit_L, N0 = sum(line$density_0),
                     K = fit_K, adjust_L = adjust_L, optim_args = optim_args)
    if (length(fit_args) > 0) for (n in names(fit_args)) fkl_args[[n]] <- fit_args[[n]]
    fkl_args[["return_optims"]] <- FALSE # << required for this fxn to work
    fits <- do.call(fit_known_leslie, fkl_args)
    if (any(is.na(fits))) {
        stop(sprintf("any(is.na(fits))\n shape = %s\n offset = %s\n K = %s\n",
                     shape, offset, K))
    }
    par_diffs <- attr(fits, "par_diffs")
    fits <- fits[param_lvls]
    par_diffs <- par_diffs[param_lvls]
    tibble(obs = c(shape, offset, K, width99(shape), med_age(shape, offset, ncol(fit_L))),
           fit = unname(fits),
           fit_diffs = unname(par_diffs),
           param = names(fits),
           rep = i) |>
        mutate(param = factor(param, levels = param_lvls))
}


# Test all `n_sims` simulation fits for a single fecundity and survival combination
test_all_sims <- function(n_sims,
                          line = line_s,
                          sigma_s = 0,
                          sigma_f = 0,
                          no_K = FALSE,
                          known_K = TRUE,
                          adjust_L = 1L,
                          fit_max_t = 30L,
                          optim_args = list(),
                          fit_args = list(max_shape = 100),
                          sim_args = list(),
                          parallel = TRUE,
                          show_progress = TRUE) {

    one_sim <- function(i) {
        fit_L <- line$leslie[,,1] |>
            random_survs(sigma_s) |>
            random_fecunds(sigma_f)
        ifits <- test_sim(i = i, line = line, fit_L = fit_L, no_K = no_K,
                          adjust_L = adjust_L, fit_max_t = fit_max_t,
                          optim_args = optim_args, fit_args = fit_args,
                          sim_args = sim_args)
        p()
        return(ifits)
    }
    p <- function() invisible(NULL)
    if (parallel) {
        if (show_progress) p <- progressor(steps = n_sims)
        out <- future_lapply(1:n_sims, one_sim, future.seed = TRUE,
                             future.packages = c("tidyverse", "aphidsync",
                                                 "gameofclones"))
    } else {
        out <- lapply(1:n_sims, one_sim)
    }
    out <- out |>
        do.call(what = bind_rows) |>
        mutate(sigma_s = sigma_s, sigma_f = sigma_f)
    return(out)
}




# ============================================================================*
# simulations ----
# ============================================================================*


overwrite <- FALSE



# Vary fecundities and survivals:
# Takes ~2 hours
if (overwrite || ! file.exists(rds_files$vary_fs)) {
    t0 <- Sys.time()
    set.seed(1013619437)
    vary_fs_fit_df <- crossing(sigma_f = c(0, sigma_f_vals),
                               sigma_s = c(0, sigma_s_vals),
                               adjust_L = 0:1) |>
        # Running it this way allows me to have an overall progress bar, not
        # one for each run:
        (\(opt_df) {
            p <- progressor(steps = nrow(opt_df))
            all_out <- opt_df |>
                pmap_dfr(\(sigma_f, sigma_s, adjust_L) {
                out <- test_all_sims(60L, sigma_s = sigma_s,
                                     sigma_f = sigma_f,
                                     known_K = known_K,
                                     adjust_L = adjust_L,
                                     optim_args = optim_args__,
                                     show_progress = FALSE) |>
                    mutate(adjust_L = adjust_L) |>
                    select(sigma_f, sigma_s, adjust_L, everything())
                p()
                return(out)
            })
            return(all_out)
        })()
    t1 <- Sys.time()
    write_rds(vary_fs_fit_df, rds_files$vary_fs)
    print(t1 - t0); rm(t0, t1)
} else vary_fs_fit_df <- read_rds(rds_files$vary_fs)



# Same thing, but just with the 2-parameter scaling
# Takes ~1 hour
if (overwrite || ! file.exists(rds_files$vary_fs2)) {
    t0 <- Sys.time()
    set.seed(966682011)
    vary_fs_fit2_df <- crossing(sigma_f = c(0, sigma_f_vals),
                                sigma_s = c(0, sigma_s_vals),
                                adjust_L = 2L) |>
        # Running it this way allows me to have an overall progress bar, not
        # one for each run:
        (\(opt_df) {
            p <- progressor(steps = nrow(opt_df))
            all_out <- opt_df |>
                pmap_dfr(\(sigma_f, sigma_s, adjust_L) {
                    out <- test_all_sims(60L, sigma_s = sigma_s,
                                         sigma_f = sigma_f,
                                         known_K = known_K,
                                         adjust_L = adjust_L,
                                         optim_args = optim_args__,
                                         show_progress = FALSE) |>
                        mutate(adjust_L = adjust_L) |>
                        select(sigma_f, sigma_s, adjust_L, everything())
                    p()
                    return(out)
                })
            return(all_out)
        })()
    t1 <- Sys.time()
    write_rds(vary_fs_fit2_df, rds_files$vary_fs2)
    print(t1 - t0); rm(t0, t1)
} else vary_fs_fit2_df <- read_rds(rds_files$vary_fs2)


vary_fs_fit_df <- bind_rows(vary_fs_fit_df, vary_fs_fit2_df)
rm(vary_fs_fit2_df)





# plots ----

vary_fs_plotter <- function(.p, .df) {

    # .df = vary_fs_fit_df; .p = "offset"; .kk = TRUE

    stopifnot(.p %in% param_lvls)
    stopifnot(.p %in% .df$param)

    pm <- list(shape = "Initial Beta shape",
               offset = "Initial Beta offset",
               K = "Density dependence",
               width99 = "Width of 99% quantile",
               med_age = "Median starting age")
    .m <- min(.df$obs, .df$fit)
    p <- .df |>
        filter(param == .p) |>
        mutate(sigma_f = factor(sigma_f, levels = sort(unique(sigma_f)),
                                labels = paste0("&sigma;<sub>f</sub> = ",
                                                sort(unique(sigma_f)))),
               sigma_s = factor(sigma_s, levels = sort(unique(sigma_s)),
                                labels = paste0("&sigma;<sub>s</sub> = ",
                                                sort(unique(sigma_s)))),
               adjust_L = factor(adjust_L, levels = 0:2,
                                 labels = paste0(0:2, "-params"))) |>
        ggplot(aes(obs, fit, color = adjust_L)) +
        # geom_hline(yintercept = .m, color = "gray70", linewidth = 0.5) +
        # geom_vline(xintercept = .m, color = "gray70", linewidth = 0.5) +
        geom_abline(slope = 1, intercept = 0, linetype = 2, color = "red") +
        geom_point() +
        xlab("Actual value") +
        ylab("Fitted value") +
        facet_wrap(~ sigma_f + sigma_s, nrow = 4, scales = "free") +
        scale_color_viridis_d(NULL, begin = 0.2, end = 0.8) +
        ggtitle(pm[[.p]]) +
        theme(strip.text = element_markdown())
    return(p)
}


vary_fs_fit_ps <- param_lvls |>
    discard(\(x) x == "K") |>
    set_names() |>
    map(vary_fs_plotter, .df = vary_fs_fit_df)


vary_fs_fit_ps[["shape"]]
vary_fs_fit_ps[["offset"]]
vary_fs_fit_ps[["width99"]]
vary_fs_fit_ps[["med_age"]]


.p <- "width99"

vary_fs_fit_df |>
    # filter(adjust_L != 2) |>
    filter(param == .p) |>
    mutate(sigma_f = factor(sigma_f, levels = sort(unique(sigma_f)),
                            labels = paste0("&sigma;<sub>f</sub> = ",
                                            sort(unique(sigma_f)))),
           sigma_s = factor(sigma_s, levels = sort(unique(sigma_s)),
                            labels = paste0("&sigma;<sub>s</sub> = ",
                                            sort(unique(sigma_s)))),
           adjust_L = factor(adjust_L, levels = 0:2,
                             labels = paste0(0:2, "-params"))) |>
    group_by(sigma_f, sigma_s, adjust_L, param) |>
    summarize(r = cor(obs, fit), .groups = "drop") |>
    ggplot(aes(r, adjust_L, color = adjust_L)) +
    geom_point(position = position_dodge(0.2)) +
    facet_wrap(~ sigma_f + sigma_s, nrow = 4, scales = "free") +
    scale_color_viridis_d(NULL, begin = 0.2, end = 0.8) +
    theme(strip.text = element_markdown())



# SUMMARY ----
#'
#' The 0-params seems to out-perform the 1- and 2-param scaling methods for
#' some reason.
#' None of them are terrible, but it's strange that the no-scaling method
#' would ever work better when the Leslie matrix isn't perfect.
#'




# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# BELOW IS NOT UPDATED ----

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



# Vary fecundities and survivals with known K:
# Takes 19.2 min
if (overwrite || ! file.exists(rds_files$vary_fs_knownK)) {
    t0 <- Sys.time()
    set.seed(1013619437)
    vary_fs_knownK_fit_df <- crossing(.surv_x = c(-1, 0, 1),
                                      .fecund_x = c(0.8, 1, 1.2)) |>
        pmap_dfr(\(.surv_x, .fecund_x) {
            out <- test_all_sims(100L, .surv_x, .fecund_x,
                                 known_K = TRUE,
                                 optim_args = optim_args__)
            cat(sprintf("finished surv_x = %.1f, fecund_x = %.1f\n",
                        .surv_x, .fecund_x))
            return(out)
        })
    t1 <- Sys.time()
    write_rds(vary_fs_knownK_fit_df, rds_files$vary_fs_knownK)
    print(t1 - t0); rm(t0, t1)
} else vary_fs_knownK_fit_df <- read_rds(rds_files$vary_fs_knownK)




# Vary fecundities and survivals, and sample every 3-4 days:
# Takes 45 min
if (overwrite || ! file.exists(rds_files$vary_fs_samp34)) {
    t0 <- Sys.time()
    set.seed(274085304)
    vary_fs_samp34_fit_df <- crossing(.surv_x = c(-1, 0, 1),
                                      .fecund_x = c(0.8, 1, 1.2)) |>
        pmap_dfr(\(.surv_x, .fecund_x) {
            out <- test_all_sims(100L, .surv_x, .fecund_x,
                                 sample_filter = TRUE,
                                 optim_args = optim_args__)
            cat(sprintf("finished surv_x = %.1f, fecund_x = %.1f\n",
                        .surv_x, .fecund_x))
            return(out)
        })
    t1 <- Sys.time()
    write_rds(vary_fs_samp34_fit_df, rds_files$vary_fs_samp34)
    print(t1 - t0); rm(t0, t1)
} else vary_fs_samp34_fit_df <- read_rds(rds_files$vary_fs_samp34)



# (old) LEFT OFF ----
#'
#' I'm trying to simulate and test fitting when reproduction is terminal only
#' (i.e., only happens in the last stage).
#' When I do this using 29 stages and repro only in the last stage, I get
#' numerical issues when I use a reproduction value in that last cell
#' that results in a leslie matrix with a lambda value equal to that from
#' the default leslie matrix.
#' For testing this idea, I therefore think I should use only 5 stages,
#' but this will require a lot of re-coding.
#'

.surv_x = -1
.fecund_x = 0.8

shape = 5.79532258792946
offset = 0.502491704541788
K = 2555.58761944115
sim_df <- sim_aphids(shape, offset, K, line_s, sd_error = FALSE,
                     sample_filter = FALSE)
sim_df

# L <- line_s$leslie[,,1]
# L[1,29] <- 10e3
# L[1,9:28] <- 0

L <- matrix(0, 5, 5)
L[1,5] <- 100
L[row(L) - col(L) == 1] <- 0.95 + 0:3*0.015
gameofclones:::sad_leslie(L)
max(abs(eigen(L, only.values = TRUE)[["values"]]))

# if (.surv_x != 0 || .fecund_x != 1) {
#     survs <- inv_logit(logit(L[row(L) - col(L) == 1]) + .surv_x)
#     L[row(L) - col(L) == 1] <- survs
#     L[1,9:29] <- L[1,9:29] * .fecund_x
# }

X <- beta_starts(shape, offset, sum(line_s$density_0), nrow(line_s$density_0))

N <- numeric(251L)
N[1] <- sum(X)

for (t in 1:250) {
    S <- 1 / (1 + N[t] / K)
    X <- S * (L %*% X)
    N[t+1] <- sum(X)
}



.line <- line_s
.line$density_0[,1] <- beta_starts(shape, offset, sum(.line$density_0),
                                   nrow(.line$density_0))
.line$density_0[,2] <- 0.0

sim_df <- gameofclones:::sim_gameofclones_full(.line,
                                               n_fields = 1,
                                               n_plants = 1,
                                               max_t = 250,
                                               plant_check_gaps = 1,
                                               wilted_N = 1e6,
                                               mean_K = K,  # 1800,
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

fits <- fit_sims(sim_df = sim_df, L = .L, K = .K,
                 compare_N = .compare_N,
                 fit_max_t = .fit_max_t, optim_args = .optim_args)






# Vary fecundities and survivals, and have only terminal reproduction:
# Takes XX min
if (overwrite || ! file.exists(rds_files$vary_fs_term)) {
    # Change `line_s` to have only terminal repro:
    line_s$leslie[1,29,1] <- 389478 # results in the same lambda
    line_s$leslie[1,9:28,1] <- 0
    t0 <- Sys.time()
    set.seed(906618959)
    vary_fs_term_fit_df <- crossing(.surv_x = c(-1, 0, 1),
                                    .fecund_x = c(0.8, 1, 1.2)) |>
        pmap_dfr(\(.surv_x, .fecund_x) {
            out <- test_all_sims(100L, .surv_x, .fecund_x,
                                 optim_args = optim_args__)
            cat(sprintf("finished surv_x = %.1f, fecund_x = %.1f\n",
                        .surv_x, .fecund_x))
            return(out)
        })
    t1 <- Sys.time()
    write_rds(vary_fs_term_fit_df, rds_files$vary_fs_term)
    print(t1 - t0); rm(t0, t1)
    # Change `line_s` back:
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
} else vary_fs_term_fit_df <- read_rds(rds_files$vary_fs_term)



#### >>>> More evaluations per box wasn't useful!! <<<<
# # Same as above but using 1000 evaluations per box:
# # Takes ~5.5 hours
# if (overwrite || ! file.exists(rds_files$vary_fs_x)) {
#     t0 <- Sys.time()
#     set.seed(1278546031)
#     vary_fs_x_fit_df <- crossing(.surv_x = c(-1, 0, 1),
#                                  .fecund_x = c(0.8, 1, 1.2)) |>
#         pmap_dfr(\(.surv_x, .fecund_x) {
#             out <- test_all_sims(100L, .surv_x, .fecund_x,
#                                  optim_args = list(box_optim = nloptr::bobyqa,
#                                                    fine_optim = optim,
#                                                    polished_optim = optim,
#                                                    n_bevals = 1000L))
#             cat(sprintf("finished surv_x = %.1f, fecund_x = %.1f\n",
#                         .surv_x, .fecund_x))
#             return(out)
#         })
#     t1 <- Sys.time()
#     write_rds(vary_fs_x_fit_df, rds_files$vary_fs_x)
#     print(t1 - t0); rm(t0, t1)
# } else vary_fs_x_fit_df <- read_rds(rds_files$vary_fs_x)
#
#
#
#
# # Same as above but comparing abundances, not growth rates
# # Takes ~49 min
# if (overwrite || ! file.exists(rds_files$vary_fs_N)) {
#     t0 <- Sys.time()
#     set.seed(106948406)
#     vary_fs_N_fit_df <- crossing(.surv_x = c(-1, 0, 1),
#                                  .fecund_x = c(0.8, 1, 1.2)) |>
#         pmap_dfr(\(.surv_x, .fecund_x) {
#             out <- test_all_sims(100L, .surv_x, .fecund_x,
#                                  compare_N = TRUE,
#                                  optim_args = optim_args__)
#             cat(sprintf("finished surv_x = %.1f, fecund_x = %.1f\n",
#                         .surv_x, .fecund_x))
#             return(out)
#         })
#     t1 <- Sys.time()
#     write_rds(vary_fs_N_fit_df, rds_files$vary_fs_N)
#     print(t1 - t0); rm(t0, t1)
# } else vary_fs_N_fit_df <- read_rds(rds_files$vary_fs_N)











# plots ----

known_fit_df |>
    split(~ rep) |>
    map_dfr(\(x) {
        shape <- c(x$obs[x$param == "shape"], x$fit[x$param == "shape"])
        offset <- c(x$obs[x$param == "offset"], x$fit[x$param == "offset"])
        m_age <- med_age(shape, offset)
        w99 <- width99(shape, offset)
        x |>
            add_row(obs = c(m_age[1], w99[1]), fit = c(m_age[2], w99[2]),
                    param = c("med_age", "width99"), rep = x$rep[1])
    }) |>
    mutate(param = factor(param, levels = param_lvls)) |>
    ggplot(aes(obs, fit)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "red") +
    geom_point() +
    facet_wrap(~ param, nrow = 1, scales = "free") +
    ggtitle("with perfect Leslie") +
    scale_color_manual(values = c("black", "dodgerblue"), guide = "none")


vary_fs_plotter <- function(.p, .df) {

    stopifnot(.p %in% param_lvls)

    if (.p %in% unique(.df$param)) {
        d <- .df |>
            filter(param == .p)
    } else {
        stopifnot(.p %in% c("width99", "p_repro", "med_age"))
        f <- eval(parse(text = .p))
        d <- .df |>
            filter(param %in% c("shape", "offset")) |>
            pivot_wider(names_from = "param", values_from = c(obs, fit)) |>
            mutate(obs = f(obs_shape, obs_offset),
                   fit = f(fit_shape, fit_offset))
    }
    g <- function(x, n) factor(x, levels = sort(unique(x)),
                               labels = paste(c("low", "perfect", "high"), n))
    pm <- list(shape = "Initial Beta shape",
               offset = "Initial Beta offset",
               K = "Density dependence",
               width99 = "Width of 99% quantile",
               med_age = "Median starting age",
               p_repro = "Proportion of initial aphids reproducing")
    .m <- min(d$obs, d$fit)
    p <- d |>
        mutate(surv_x = g(surv_x, "survival"),
               fecund_x = g(fecund_x, "fecundity")) |>
        ggplot(aes(obs, fit)) +
        # geom_hline(yintercept = .m, color = "gray70", linewidth = 0.5) +
        # geom_vline(xintercept = .m, color = "gray70", linewidth = 0.5) +
        geom_abline(slope = 1, intercept = 0, linetype = 2, color = "red") +
        geom_point() +
        xlab("Actual value") +
        ylab("Fitted value") +
        facet_wrap(~ surv_x + fecund_x, nrow = 3, scales = "free") +
        ggtitle(pm[[.p]])
    return(p)
}


vary_fs_fit_ps <- param_lvls |>
    set_names() |>
    map(vary_fs_plotter, .df = vary_fs_fit_df)

# vary_fs_fit_ps[["shape"]]
vary_fs_fit_ps[["width99"]]
vary_fs_fit_ps[["med_age"]]
vary_fs_fit_ps[["offset"]]
vary_fs_fit_ps[["K"]]

vary_fs_samp34_fit_ps <- param_lvls |>
    set_names() |>
    map(vary_fs_plotter, .df = vary_fs_samp34_fit_df)

vary_fs_samp34_fit_ps[["width99"]] +
vary_fs_samp34_fit_ps[["med_age"]] +
vary_fs_samp34_fit_ps[["offset"]] +
vary_fs_samp34_fit_ps[["K"]] +
    plot_layout(nrow = 2) &
    theme(strip.text = element_text(size = 6),
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 6))


# for (p in c("width99", "offset", "K", "p_repro")) {
#     ggsave(sprintf("_figures/fit-%s.pdf", p), vary_fs_fit_ps[[p]],
#            width = 6, height = 6)
# }; rm(p)





vary_fs_N_fit_ps <- param_lvls |>
    set_names() |>
    map(vary_fs_plotter, .df = vary_fs_N_fit_df)

# vary_fs_N_fit_ps[["shape"]]
vary_fs_N_fit_ps[["width99"]]
vary_fs_N_fit_ps[["offset"]]
vary_fs_N_fit_ps[["K"]]
vary_fs_N_fit_ps[["p_repro"]]


# for (p in c("width99", "offset", "K", "p_repro")) {
#     ggsave(sprintf("_figures/fit-N-%s.pdf", p), vary_fs_N_fit_ps[[p]],
#            width = 6, height = 6)
# }; rm(p)




vary_fs_knownK_fit_ps <- param_lvls |>
    set_names() |>
    map(vary_fs_plotter, .df = vary_fs_knownK_fit_df)

# vary_fs_knownK_fit_ps[["shape"]]
vary_fs_knownK_fit_ps[["width99"]]
vary_fs_fit_ps[["med_age"]]
vary_fs_knownK_fit_ps[["offset"]]
vary_fs_knownK_fit_ps[["K"]]



# for (p in c("width99", "offset", "K", "p_repro")) {
#     ggsave(sprintf("_figures/fit-X-%s.pdf", p), vary_fs_x_fit_ps[[p]],
#            width = 6, height = 6)
# }; rm(p)






