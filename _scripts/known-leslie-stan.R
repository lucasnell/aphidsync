

source("_scripts/known-leslie-simple-preamble.R")

suppressPackageStartupMessages({
    library(gameofclones)
    library(rstan)
})



options(mc.cores = max(1L, parallel::detectCores() - 2L))
rstan_options(auto_write = TRUE)
rstan_options(threads_per_chain = 2)


no_scaling_mod <- stan_model(file = "_scripts/no-scaling.stan",
                             model_name = "no-scaling")

with_scaling_mod <- stan_model(file = "_scripts/1par-scaling.stan",
                               model_name = "1par-scaling")


n_stages = 10L
n_repros = 5L
sigma_s = sigma_s_vals[[2]]
sigma_f = sigma_f_vals[[2]]


set.seed(1073136633); {
    sim_L <- make_rand_L(n_stages, n_repros)
    fit_L <- sim_L |>
        random_survs(sigma_s) |>
        random_fecunds(sigma_f)
    fit_L[1,] <- fit_L[1,] * 0.9
    obs <- list(shape = runif(1, 1, 6),
                offset = runif(1, 0, 1),
                sigma = 0.1)
    obs$med_age <- med_age(obs$shape, obs$offset, n_stages)
    obs$width99 = width99(obs$shape)
    K <- runif(1, 1000, 5000)
    re_df <- sim_pop(obs$shape, obs$offset, K, sim_L) |>
        arrange(time) |>
        mutate(re = (lead(N) / N)^(1/(lead(time) - time)),
               N2 = N * exp(rnorm(n(), 0, obs$sigma)),
               re2 = (lead(N2) / N2)^(1/(lead(time) - time)),
               re3 = re + rnorm(n(), 0, obs$sigma)) |>
        head(30L)
    data_list <- list(n_obs = nrow(re_df),
                      n_stages = n_stages,
                      total_N0 = re_df$N[1],
                      K = K,
                      L = fit_L,
                      pcg = re_df$re3,
                      time = re_df$time + 1L,
                      shape_mean = 3.5,
                      shape_sd = 3,
                      offset_mean = 0.5,
                      offset_sd = 1/12,
                      sigma_mean = 0.1,
                      sigma_sd = 0.1)
}


fit0 <- sampling(no_scaling_mod,
                 data = data_list,
                 chains = 4, iter = 4000, control = list(max_treedepth = 20))
fit1 <- sampling(with_scaling_mod,
                 data = c(data_list, list(fx_sd = 1)),
                 chains = 4, iter = 4000, control = list(max_treedepth = 20))


int_pars <- c(param_lvls[param_lvls != "K"], "sigma")


# pairs(fit0, pars = int_pars)
# pairs(fit1, pars = int_pars)
# pairs(fit2, pars = int_pars)

plot(fit0, plotfun = "trace", pars = int_pars)
plot(fit1, plotfun = "trace", pars = int_pars)

real_compare_plotter <- function(fit, .title = waiver()) {
    plot(fit, plotfun = "hist", pars = int_pars, bins = 30) +
        geom_vline(data = tibble(parameter = factor(int_pars),
                                 vline = obs[int_pars] |> unlist()),
                   aes(xintercept = vline), linetype = 1, color = "dodgerblue",
                   linewidth = 2) +
        labs(title = .title)
}

real_compare_plotter(fit0, "no scaling")
real_compare_plotter(fit1, "1-par scaling")

cbind(obs = unlist(obs[int_pars]),
      no_scale = get_posterior_mean(fit0, int_pars)[,"mean-all chains"],
      scale = get_posterior_mean(fit1, int_pars)[,"mean-all chains"])


tibble(time = re_df$time,
       obs = re_df$re,
       no_scale = unname(get_posterior_mean(fit0, "pcg_pred")[,"mean-all chains"]),
       scale = unname(get_posterior_mean(fit1, "pcg_pred")[,"mean-all chains"])) |>
    (\(x) {
        cat("cor(obs, no_scale) = ", cor(x$obs, x$no_scale), "\n")
        cat("cor(obs, scale) = ", cor(x$obs, x$scale), "\n")
        return(x)
    })() |>
    pivot_longer(no_scale:scale, names_to = "method", values_to = "fit") |>
    ggplot(aes(time, obs)) +
    geom_line(aes(y = fit), color = "gray70", linewidth = 1) +
    geom_line(color = "black", linetype = "22", linewidth = 0.75) +
    facet_wrap(~ method, nrow = 1)



