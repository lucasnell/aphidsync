
suppressPackageStartupMessages({
    library(tidyverse)
    library(aphidsync)
    library(future.apply)
    library(progressr)
    library(patchwork)
    library(ggtext)
    library(rstan)
})

plan(multisession, workers = max(parallel::detectCores() - 2L, 1L))
options(mc.cores = max(1L, parallel::detectCores() - 2L))
rstan_options(auto_write = TRUE)
rstan_options(threads_per_chain = 2)
handlers(global = TRUE)
handlers("progress")


if (interactive() && file.exists(".Rprofile")) source(".Rprofile")




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

# List of Leslie matrices by number of stages where repro occurs:
LL <- rep(list(NA), 4)
LL[[1]] <- matrix(0, 5, 5)
LL[[1]][row(LL[[1]]) - col(LL[[1]]) == 1] <- 0.96 * exp(0.2 * 0:-3)
LL[[1]][1,5] <- 10

for (i in 2:length(LL)) {
    LL[[i]] <- LL[[1]]
    f_inds <- (nrow(LL[[1]]) - i + 1L):nrow(LL[[1]])
    # Make fecundities decrease with age:
    LL[[i]][1,f_inds] <- length(f_inds):1
    # Use optimize to scale these values so that the growth rate of this
    # Leslie matrix closely matches that for `LL[[1]]`:
    LL[[i]][1,f_inds] <- LL[[i]][1,f_inds] * optimize(\(x) {
        fL <- LL[[i]]
        fL[1,f_inds] <- fL[1,f_inds] * x
        abs(calc_lambda(fL) - calc_lambda(LL[[1]]))
    }, c(0, LL[[1]][1,5]), tol = .Machine$double.eps^0.5)[["minimum"]]
}

sapply(LL, calc_lambda)




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




# =============================================================================*
# fit testing ----
# =============================================================================*


idx <- 4L
sim_L <- LL[[idx]]
fit_L <- adjust_survs(LL[[idx]], 0.2)|>
    adjust_fecunds(0.1)

set.seed(78); {
    shape <- runif(1, 1, 6)
    offset <- runif(1, 0, 1)
    sigma <- 0.1
    ma <- med_age(shape, offset, 5L)
    w99 <- width99(shape)
    K <- runif(1, 1000, 5000)
    re_df <- sim_aphids(shape, offset, K, sim_L) |>
        arrange(time) |>
        mutate(re = (lead(N) / N)^(1/(lead(time) - time)),
               N2 = N * exp(rnorm(n(), 0, sigma)),
               re2 = (lead(N2) / N2)^(1/(lead(time) - time))) |>
        head(30L)
}


fit <- stan(file = "_scripts/no-scaling.stan",
            model_name = "no-scaling",
            data = list(n_obs = nrow(re_df),
                        n_stages = 5L,
                        total_N0 = 100,
                        K = K,
                        L = fit_L,
                        shape_mu = 3.5,
                        shape_sigma = 3,
                        pcg = re_df$re2,
                        time = re_df$time + 1L),
            chains = 4, iter = 4000, control = list(max_treedepth = 20))


int_pars <- c("shape", "offset", "med_age", "width99")

# pairs(fit, pars = int_pars)
plot(fit, plotfun = "trace", pars = int_pars)
plot(fit, plotfun = "hist", pars = int_pars, bins = 30) +
    geom_vline(data = tibble(parameter = factor(int_pars),
                             vline = c(shape, offset, ma, w99)),
               aes(xintercept = vline), linetype = 1, color = "dodgerblue",
               linewidth = 2)



