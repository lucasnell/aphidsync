
suppressPackageStartupMessages({
    library(Hmisc)  # must come first
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

param_lvls <- c("shape", "offset", "K", "width99", "med_age")



# Make list of Leslie matrices by number of stages where repro occurs:
LL <- rep(list(NA), 4)
LL[[1]] <- matrix(0, 5, 5)
LL[[1]][row(LL[[1]]) - col(LL[[1]]) == 1] <- 0.96 * exp(0.2 * 0:-3)
LL[[1]][1,5] <- 10

for (i in 2:length(LL)) {
    LL[[i]] <- LL[[1]]
    f_inds <- (nrow(LL[[1]]) - i + 1L):nrow(LL[[1]])
    # Make fecundities decrease with age:
    LL[[i]][1,f_inds] <- length(f_inds):1
    # Scale these values so that the growth rate of this
    # Leslie matrix closely matches that for `LL[[1]]`:
    LL[[i]] <- adjust_lambda(LL[[i]], calc_lambda(LL[[1]]), "f",
                             LL[[1]][1,5], direct_l = TRUE, optim_tol = 0)
}; rm(i, f_inds)

# sapply(LL, calc_lambda)




# =============================================================================*
# functions ----
# =============================================================================*

sim_pop <- function(shape, offset, K, L,
                    sample_filter = FALSE,
                    max_t = 150,
                    sd_error = 0) {

    abunds0 <- beta_starts(shape, offset, 100, nrow(L))
    .time <- 0:max_t

    sim_df <- tibble(N = as.numeric(sim_N(abunds0, L, .time, K)),
                     time = .time)

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
test_sim <- function(i, sim_L, fit_L,
                     no_K = FALSE,
                     adjust_L = 0L,
                     fit_max_t = 30L,
                     optim_args = list(),
                     fit_args = list(),
                     sim_args = list()) {
    shape <- runif(1, 1, 6)
    offset <- runif(1, 0, 1)
    if (no_K) {
        K <- Inf
    } else K <- runif(1, 1000, 5000)
    sp_args <- list("shape" = shape, "offset" = offset, "K" = K, "L" = sim_L)
    if (length(sim_args) > 0) for (n in names(sim_args)) sp_args[[n]] <- sim_args[[n]]
    sim_df <- do.call(sim_pop, sp_args)
    fkl_args <- list(sim_df = sim_df, L = fit_L, N0 = 100, K = K,
                     adjust_L = adjust_L, optim_args = optim_args)
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
    tibble(obs = c(shape, offset, K, width99(shape), med_age(shape, offset, ncol(sim_L))),
           fit = unname(fits),
           fit_diffs = unname(par_diffs),
           param = names(fits),
           rep = i) |>
        mutate(param = factor(param, levels = param_lvls))
}


# Make a random Leslie matrix:
make_rand_L <- function(n_stages, n_repros) {
    # -----------*
    # Randomly generated items
    # -----------*
    # Shape and scale for Weibull distribution used to generate
    # fecundities (Weibull used bc it fit best for aphid data)
    # Both distributions are uniforms with means that are about what the aphid
    # Leslie matrix is
    weib_shape <- runif(1, 0, 4)
    weib_scale <- runif(1, 0, 16)
    # Survivals:
    survs <- inv_logit(rnorm(n_stages-1, logit(0.9), 1))
    # Overall growth rate:
    lambda <- runif(1, 1.05, 1.2)
    # -----------*
    # Creating Leslie matrix
    # -----------*
    L <- make_L1(weib_shape, weib_scale, n_stages, n_stages - n_repros + 1L)
    L[row(L) - col(L) == 1] <- survs
    L <- adjust_lambda(L, lambda)
    return(L)
}


# Test all `n_sims` simulation fits for a single fecundity and survival combination
# If L is NULL and both n_stages and n_repros are provided, it creates a new,
# randomly generated leslie matrix for each iteration.
test_all_sims <- function(n_sims, L,
                          sigma_s = 0,
                          sigma_f = 0,
                          n_stages = NULL,
                          n_repros = NULL,
                          no_K = FALSE,
                          adjust_L = 0L,
                          fit_max_t = 30L,
                          optim_args = list(),
                          fit_args = list(max_shape = 100),
                          sim_args = list(),
                          parallel = TRUE,
                          show_progress = TRUE) {
    if (is.null(L) && (is.null(n_stages) || is.null(n_repros))) {
        stop("If L is NULL, n_stages and n_repros must both be provided.")
    }
    one_sim <- function(i) {
        if (is.null(L)) {
            sim_L <- make_rand_L(n_stages, n_repros)
        } else {
            sim_L <- L
        }
        fit_L <- sim_L |>
            random_survs(sigma_s) |>
            random_fecunds(sigma_f)
        ifits <- test_sim(i = i, sim_L = sim_L, fit_L = fit_L, no_K = no_K,
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
                             future.packages = c("tidyverse", "aphidsync"))
    } else {
        out <- lapply(1:n_sims, one_sim)
    }
    out <- out |>
        do.call(what = bind_rows) |>
        mutate(sigma_s = sigma_s, sigma_f = sigma_f)
    return(out)
}



# These two combinations look good for testing:
sigma_s_vals <- c(0.1, 0.2, 0.5)
sigma_f_vals <- c(0.05, 0.1, 0.2)

# compare_f <- function(nsims, sigma_f) {
#     L <- matrix(0, 5, 5)
#     out <- map_dfr(1:nsims, \(i) {
#         L[1,] <- runif(5, 0, 10)
#         L2 <- adjust_fecunds(L, sigma_f)
#         return(tibble(rep = i, sigma_f = sigma_f,
#                       before = L[1,], after = L2[1,]))
#     })
#     return(out)
# }
# compare_s <- function(nsims, sigma_s) {
#     L <- matrix(0, 5, 5)
#     s_inds <- which(row(L) - col(L) == 1)
#     out <- map_dfr(1:nsims, \(i) {
#         L[s_inds] <- runif(length(s_inds))
#         L2 <- adjust_survs(L, sigma_s)
#         return(tibble(rep = i, sigma_s = sigma_s,
#                       before = L[s_inds], after = L2[s_inds]))
#     })
#     return(out)
# }
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


