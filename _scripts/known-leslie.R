
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


# Files created here to store simulation data:
rds_files <- list(known = "_data/known_fit_df.rds",
                  vary_fs = "_data/vary_fs_fit_df.rds",
                  vary_fs_knownK = "_data/vary_fs_knownK_fit_df.rds",
                  vary_fs_N = "_data/vary_fs_N_fit_df.rds")


col_pal <- c(resistant = viridis(100)[50],
             susceptible = viridis(100)[95],
             wasp = viridis(100)[1])

param_lvls <- c("shape", "offset", "K", "med_age", "width99")

# Shared options for `winnowing_optim`:
optim_args__ <- list(box_optim = nloptr::bobyqa,
                     fine_optim = nloptr::bobyqa,
                     polished_optim = nloptr::bobyqa,
                     box_control = list(maxit = 500L, reltol = 1e-6),
                     fine_control = list(maxit = 1000L, reltol = 1e-8),
                     polished_control = list(maxit = 10e3L, reltol = 1e-10),
                     n_fine = 500L,
                     n_polished = 250L)

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





sim_aphids <- function(.shape, .offset, .K, .line,
                       sample_filter = FALSE,
                       sd_error = 0,
                       .max_t = 250) {

    .line$density_0[,1] <- beta_starts(.shape, .offset, sum(.line$density_0),
                                       nrow(.line$density_0))
    .line$density_0[,2] <- 0.0

    sim_df <- gameofclones:::sim_gameofclones_full(.line,
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

    if (sd_error > 0) {
        sim_df$N <- sim_df$N * exp(rnorm(nrow(sim_df), 0, sd_error))
    }

    return(sim_df)

}





#
# Fit the two shape parameters for a set of simulations.
# NOTE: This version doesn't allow for multiple lines (yet).
#
fit_sims <- function(sim_df, L, K, max_shape = 800,
                     compare_N = FALSE,
                     fit_max_t = 30L,
                     optim_args = list()) {
                     # cycles = TRUE,

    # sim_df = sim_aphids(2, 0, 1800, line_s); L = line_s$leslie[,,1]
    # max_shape = 800; compare_N = FALSE; optim_args = list()

    stopifnot(length(unique(sim_df$line)) == 1L)

    if (nrow(sim_df) < 3L) return("low rows") # special meaning in `test_sim`

    time <- sim_df$time[1:fit_max_t]
    if (compare_N) {
        observed <- sim_df$N[1:fit_max_t]
    } else {
        re_df <- sim_df |>
            arrange(time) |>
            mutate(re = (lead(N) / N)^(1/(lead(time) - time)))
        observed <- re_df$re[1:fit_max_t]
    }

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
    #
    # # Fit across multiple starting offsets, then choose best-fitting one:
    # val <- 1e10
    # op <- NULL
    # n_starts <- 100L
    # guess_list <- lapply(1:n_starts, \(i) runif(3, 0, c(10, 1, 10e3)))
    #
    # # Choose 10 random reps to have very low and high offset values.
    # # This helps with fitting.
    # n_hilo <- 10L
    # off_inds <- sample.int(n_starts, n_hilo)
    # for (i in off_inds[1:(n_hilo %/% 2L)])
    #     guess_list[[i]][2] <- runif(1, 0, 1e-6)
    # for (i in off_inds[(n_hilo %/% 2L + 1L):n_hilo])
    #     guess_list[[i]][2] <- runif(1, 1 - 1e-6, 1)
    # for (s in 1:n_starts) {
    #     guess <- guess_list[[s]]
    #     new_op <- optim(guess, known_fit_aphids0,
    #                     L = L, re = re_df$re, time = re_df$time,
    #                     max_shape = max_shape)
    #     if (new_op$convergence == 0 && new_op$value < val) {
    #         val <- new_op$value
    #         op <- new_op
    #     }
    #     if (new_op$convergence == 0 && !is.null(op) && new_op$value == val) {
    #         if (!isTRUE(all.equal(op$par, new_op$par))) {
    #             # This triggers a specific error in `test_sim`:
    #             return("multiple models")
    #         }
    #     }
    # }
    #
    # if (is.null(op)) return("no converged") # specific meaning in `test_sim`

    # If not provided as a known value, estimate K from time series and
    # Leslie matrix:
    if (is.na(K)) {
        pred_K <- median(sim_df$N[sim_df$time > 100]) /
            (max(abs(eigen(L)[["values"]])) - 1)
        if (is.na(pred_K)) stop("is.na(pred_K)")
    } else pred_K <- K

    # Do optimization using winnowing approach
    op_args <- optim_args
    op_args[["fn"]] <- known_fit_aphids0
    op_args[["lower_bounds"]] <- c(0,   0)
    op_args[["upper_bounds"]] <- c(100, 1)
    op_args[["fn_args"]] <- list(L = L,
                                 K = pred_K,
                                 obs = observed,
                                 time = time,
                                 max_shape = max_shape,
                                 compare_N = compare_N)
    ops <- do.call(winnowing_optim, op_args)
    op_pars <- sapply(ops, \(x) x$par)
    op_diffs <- pmin(apply(op_pars, 1, \(x) diff(range(x)) / median(x)),
                     apply(op_pars, 1, \(x) diff(range(x))))
    if (max(op_diffs) > 0.01) warning("max(op_diffs) > 0.01")
    op <- ops[[1]]

    if (!"par" %in% names(op)) stop("unknown output parameter name")
    pars <- c(op$par, pred_K)
    names(pars) <- c("shape", "offset", "K")

    return(pars)

}



# Test 1 simulation fit for a single fecundity and survival combination
test_sim <- function(i, .L, .sd_error, .known_K, .compare_N, .fit_max_t,
                     .optim_args = list()) {
    shape <- runif(1, 0.5, 6)
    offset <- runif(1, 0, 1)
    K <- runif(1, 1000, 5000)
    .K <- ifelse(.known_K, K, NA_real_)
    # Simulate `line_s`, but fit using altered `.line`:
    sim_df <- sim_aphids(shape, offset, K, line_s, sd_error = .sd_error)
    fits <- fit_sims(sim_df = sim_df, L = .L, K = .K,
                     compare_N = .compare_N,
                     fit_max_t = .fit_max_t, optim_args = .optim_args)
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
    tibble(obs = c(shape, offset, K),
           fit = unname(fits),
           param = names(fits),
           rep = i) |>
        mutate(param = factor(param, levels = param_lvls))
}


# Test all 100 simulation fits for a single fecundity and survival combination
test_all_sims <- function(n_sims, .surv_x = 0, .fecund_x = 1, sd_error = 0,
                          known_K = FALSE,
                          compare_N = FALSE, fit_max_t = 30L,
                          optim_args = list()) {
    L <- line_s$leslie[,,1]
    if (.surv_x != 0 || .fecund_x != 1) {
        survs <- inv_logit(logit(L[row(L) - col(L) == 1]) + .surv_x)
        L[row(L) - col(L) == 1] <- survs
        L[1,9:29] <- L[1,9:29] * .fecund_x
    }
    p <- progressor(steps = n_sims)
    out <- future_lapply(1:n_sims, function(i) {
        ifits <- test_sim(i, .L = L, .sd_error = sd_error, .known_K = known_K,
                          .compare_N = compare_N, .fit_max_t = fit_max_t,
                          .optim_args = optim_args)
        p()
        return(ifits)
    },
    future.seed = TRUE,
    future.packages = c("tidyverse", "gameofclones", "aphidsync")) |>
        do.call(what = bind_rows) |>
        mutate(surv_x = .surv_x, fecund_x = .fecund_x)
    return(out)
}






# simulations ----

overwrite <- FALSE

# Takes ~2.3 min
if (overwrite || ! file.exists(rds_files$known)) {
    t0 <- Sys.time()
    set.seed(1740563097)
    known_fit_df <- test_all_sims(100L, optim_args = optim_args__) |>
        select(-surv_x, -fecund_x)
    t1 <- Sys.time()
    write_rds(known_fit_df, rds_files$known)
    print(t1 - t0); rm(t0, t1)
} else known_fit_df <- read_rds(rds_files$known)





# Vary fecundities and survivals:
# Takes 19.2 min
if (overwrite || ! file.exists(rds_files$vary_fs)) {
    t0 <- Sys.time()
    set.seed(1013619437)
    vary_fs_fit_df <- crossing(.surv_x = c(-1, 0, 1),
                               .fecund_x = c(0.8, 1, 1.2)) |>
        pmap_dfr(\(.surv_x, .fecund_x) {
            out <- test_all_sims(100L, .surv_x, .fecund_x,
                                 optim_args = optim_args__)
            cat(sprintf("finished surv_x = %.1f, fecund_x = %.1f\n",
                        .surv_x, .fecund_x))
            return(out)
        })
    t1 <- Sys.time()
    write_rds(vary_fs_fit_df, rds_files$vary_fs)
    print(t1 - t0); rm(t0, t1)
} else vary_fs_fit_df <- read_rds(rds_files$vary_fs)



# LEFT OFF ----

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






