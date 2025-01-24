
suppressPackageStartupMessages({
    library(tidyverse)
    library(viridisLite)
    library(gameofclones)
    library(parallel)
})

options(mc.cores = max(detectCores() - 2L, 1L))


# ggplot2 theme:
theme_set(theme_classic() %+replace%
              theme(strip.background = element_blank(),
                    strip.text = element_text(size = 11),
                    legend.background = element_blank(),
                    plot.title = element_text(size = 14, hjust = 0)))

if (file.exists(".Rprofile")) source(".Rprofile")



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


# calculate beta mean and variance from shapes:
beta_mean <- function(shapes) {
    shapes[1] / (shapes[1] + shapes[2])
}
beta_var <- function(shapes) {
    (shapes[1] * shapes[2]) / ((shapes[1] + shapes[2])^2 *
                                   (shapes[1] + shapes[2] + 1))
}

beta_starts <- function(shape1, shape2, total_aphids0, compartments = 29L) {
    stopifnot(shape1 > 0 && shape2 > 0)
    stopifnot(compartments > 0 && total_aphids0 > 0)
    stopifnot(compartments == round(compartments))
    stopifnot(total_aphids0 == round(total_aphids0))
    times = seq(0, 1, length.out = round(compartments + 1))
    pbeta_vals = pbeta(times, shape1 = shape1, shape2 = shape2)
    aphids0 = total_aphids0 * diff(pbeta_vals)
    if(round(sum(aphids0)) != round(total_aphids0)){
        warning(paste("beta_starts magnitude error (",
                      round(sum(aphids0)) - round(total_aphids0), ", shape1 = ",
                      round(shape1, digits = 2), ", shape2 = ",
                      round(shape2, digits = 2), sep = ''), immediate. = T)}
    if(round(length(aphids0)) != round(compartments)){
        warning(paste("beta_starts length error, shape1 = ",
                      round(shape1, digits = 2), ", shape2 = ",
                      shape2, sep=''), immediate. = T)}
    return(aphids0)
}



sim_aphids <- function(.shape1, .shape2, sample_filter = FALSE) {

    line <- line_s
    line$density_0[,1] <- beta_starts(.shape1, .shape2, sum(line$density_0),
                                      nrow(line$density_0))
    line$density_0[,2] <- 0.0

    sim_df <- gameofclones:::sim_gameofclones_full(line,
                                                   n_fields = 1,
                                                   n_plants = 16,
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
                                                   mean_K = 1800,
                                                   clear_surv = 0,
                                                   wasp_density_0 = 0,
                                                   wasp_delay = 8,
                                                   mum_smooth = 0,
                                                   aphid_plant_disp_p = 0.1,
                                                   plant_disp_mort = 0,
                                                   wilted_mort = 0.3,
                                                   extinct_N = 1) |>
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
fit_aphids0 <- function(pars, Lvec, re_df, log_trans, .plot = FALSE) {

    # pars <- c(533.764, 393.973)

    if (log_trans) {
        .shape1 <- exp(pars[[1]])
        .shape2 <- exp(pars[[2]])
    } else {
        if (any(pars < 0)) return(1e10)
        .shape1 <- pars[[1]]
        .shape2 <- pars[[2]]
    }

    aphids0 <- beta_starts(.shape1, .shape2, total_aphids0 = 32)

    # Extrapolate through time:
    Nmat  <- matrix(NA_real_, nrow = length(Lvec) + 1, ncol = length(aphids0))
    Nmat[1, ] <- aphids0
    for (i in 2:(nrow(Nmat))) {
        Nmat[i, ] = Lvec[[i-1]] %*% Nmat[(i - 1), ]
    }

    sumNvec <- apply(Nmat, 1, sum)
    re_pred <- (lead(sumNvec) / sumNvec)[(re_df$time + 1L)]

    if (.plot) {
        plot(re_pred, ylim = c(0, max(c(re_pred, re_df$re), na.rm = TRUE)),
             type = "l", ylab = "Re predicted (black) and observed (red)",
             xlab = "time", lwd = 2,
             main = sprintf("shape1 = %.4f, shape2 = %.4f", .shape1, .shape2))
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
fit_sims <- function(sim_df, cycles = TRUE, log_trans = FALSE, real_L = TRUE) {

    stopifnot(length(unique(sim_df$line)) == 1L)

    if (nrow(sim_df) < 3L) return(c(NA_real_, NA_real_))

    n_stages <- formals(beta_starts)[["compartments"]]
    stopifnot(n_stages %in% c(14, 29))
    adult_stage <- ifelse(n_stages == 14, 5, 9) # first adult stage

    # Setup starting leslie matrix with only survivals = 1 (without fecundities)
    leslie <- matrix(0, n_stages, n_stages)
    # survivals (set to 1 for now):
    leslie[(row(leslie) - col(leslie) == 1)] <- 1
    # fecundities (known):
    leslie[1, adult_stage:n_stages] <- line_s$leslie[1, adult_stage:n_stages,1]

    re_df <- sim_df |>
        arrange(time) |>
        mutate(re = (lead(N) / N)^(1/(lead(time) - time)))

    amr_df <- sim_df
    if (cycles) {
        amr_df <- amr_df |>
            filter(time %% 29L == 0)
    }
    amr_df <- amr_df |>
        mutate(amr = (lead(N) / N)^(1/(lead(time) - time))) |>
        filter(!is.na(amr))


    amr_f <- approxfun(amr_df$time, y = amr_df$amr, yleft = amr_df$amr[1],
                       yright = tail(amr_df$amr, 1))

    amr_pred = amr_f(1:max(re_df$time))

    if (real_L) {
        Lvec <- lapply(amr_pred, \(amr) {
            # op <- optimize(function(x) {
            #     L <- line_s$leslie[,,1]
            #     # L[1, adult_stage:n_stages] <- L[1, adult_stage:n_stages] * x
            #     L[(row(L) - col(L) == 1)] <- L[(row(L) - col(L) == 1)] * x
            #     lam <- max(abs(eigen(L, only.values = TRUE)[["values"]]))
            #     return(abs(lam - amr))
            # }, c(0, amr), tol = 1e-8)
            # L <- line_s$leslie[,,1]
            # # L[1, adult_stage:n_stages] <- L[1, adult_stage:n_stages] *
            # #     op$minimum
            # L[(row(L) - col(L) == 1)] <- L[(row(L) - col(L) == 1)] * op$minimum
            L <- line_s$leslie[,,1]
            return(L)
        })
    } else {
        Lvec <- lapply(amr_pred, \(amr) {
            L_fit_fun <- function(x) {
                L <- leslie
                L[1, adult_stage:n_stages] <- x[1]
                L[(row(L) - col(L) == 1)] <- x[2]
                lam <- max(abs(eigen(L, only.values = TRUE)[["values"]]))
                return(abs(lam - amr))
            }
            status <- 1
            iters <- 0
            while (status != 0 && iters < 10) {
                op <- optim(c(1, 1), L_fit_fun)
                status <- op$convergence
                first_guess <- runif(2, 0, 10)
                iters <- iters + 1
            }
            if (status != 0) stop("Leslie could not converge")
            L <- leslie
            L[1, adult_stage:n_stages] <- op$par[1]
            L[(row(L) - col(L) == 1)] <- op$par[2]
            return(L)
        })
    }

    if (log_trans) {
        first_guess <- c(0, 0)
    } else first_guess <- c(1, 1)
    status <- 1
    iters <- 0
    while (status != 0 && iters < 10) {
        op <- optim(first_guess, fit_aphids0,
                    Lvec = Lvec, re_df = re_df, log_trans = log_trans)
        status <- op$convergence
        first_guess <- runif(2, 0, 10)
        iters <- iters + 1
    }

    if (status != 0) return(c(NA_real_, NA_real_))

    if (log_trans) {
        pars <- exp(op$par)
    } else pars <- op$par

    return(pars)

}



# LEFT OFF -----
# not working very well...

# fit_aphids0(c(2.394014, 5.196421))
# fit_aphids0(c(2.394014, 5.196421), TRUE)
# fit_aphids0(c(2, 5), TRUE)

# Takes ~1 min total
{
    t0 <- Sys.time()
    RNGkind("L'Ecuyer-CMRG")
    set.seed(1740563097)
    fit_df <- mclapply(1:100, \(x) {
        shape1 <- runif(1, 0.5, 6)
        shape2 <- runif(1, 0.5, 9)
        obs_mean <- beta_mean(c(shape1, shape2))
        obs_var <- beta_var(c(shape1, shape2))
        sim_df <- sim_aphids(shape1, shape2)
        fits <- fit_sims(sim_df, real_L = TRUE)
        fit_mean <- beta_mean(fits)
        fit_var <- beta_var(fits)
        tibble(obs = c(shape1, shape2, obs_mean, obs_var),
               fit = c(fits, fit_mean, fit_var),
               param = c("shape1", "shape2", "mean", "var")) |>
            mutate(param = factor(param, levels = param))
    }) |>
        imap_dfr(\(x, i) {
            mutate(x, rep = i)
        })
    cat("Finished #1!\n")
    t1 <- Sys.time()
    print(t1 - t0)
    t0 <- Sys.time()
    set.seed(1376453977)
    fit_df2 <- mclapply(1:100, \(x) {
        shape1 <- runif(1, 0.5, 6)
        shape2 <- runif(1, 0.5, 9)
        obs_mean <- beta_mean(c(shape1, shape2))
        obs_var <- beta_var(c(shape1, shape2))
        sim_df <- sim_aphids(shape1, shape2)
        fits <- fit_sims(sim_df, real_L = FALSE)
        fit_mean <- beta_mean(fits)
        fit_var <- beta_var(fits)
        tibble(obs = c(shape1, shape2, obs_mean, obs_var),
               fit = c(fits, fit_mean, fit_var),
               param = c("shape1", "shape2", "mean", "var")) |>
            mutate(param = factor(param, levels = param))
    }) |>
        imap_dfr(\(x, i) {
            mutate(x, rep = i)
        })
    cat("Finished #2!\n")
    print(t1 - t0); rm(t0, t1)
    RNGkind("default")
}




fit_df |>
    # filter(fit < 15) |>
    ggplot(aes(obs, fit)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "red") +
    geom_point() +
    facet_wrap(~ param, nrow = 2, scales = "free") +
    ggtitle("using known Leslie matrix") +
    # coord_equal() +
    scale_y_continuous()

fit_df2 |>
    # filter(fit < 15) |>
    ggplot(aes(obs, fit)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "red") +
    geom_point() +
    facet_wrap(~ param, nrow = 2, scales = "free") +
    ggtitle("using unknown Leslie matrix") +
    # coord_equal() +
    scale_y_continuous()


