
#'
#' This script uses a polynomial regression to predict what to multiply
#' fecundity values in a Leslie matrix by to get a lambda = 1.291709.
#' This lambda is the value of the Leslie matrix used in the simulations.
#'
#' These fecundity values are fit to a Weibull distribution of a given shape
#' and scale and are standardized to sum to 1.
#'.
#'


suppressPackageStartupMessages({
    library(tidyverse)
    library(viridisLite)
    library(gameofclones)
    library(aphidsync)
    library(future.apply)
    library(progressr)
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

lambda__ <- max(abs(eigen(line_s$leslie[,,1], only.values = TRUE)[["values"]]))


# # This will inform `get_x` below to see what the upper bound should be
# # (it's costly to just set it to a very high number for all reps):
# # Takes ~3 sec
# get_min_x <- function(.shape, .scale) {
#     stopifnot(length(.shape) == length(.scale))
#     stopifnot(all(.shape > 0))
#     stopifnot(all(.scale > 0))
#     n_sims <- length(.shape)
#     p <- progressor(steps = n_sims)
#     out <- future_mapply(\(sh, sc) {
#         L <- matrix(0, 29, 29)
#         L[(row(L) - col(L) == 1)] <- 1
#         p_distr_vals <- pweibull(0:21, shape = sh, scale = sc)
#         L[1, 9:29] <- diff(p_distr_vals) / sum(diff(p_distr_vals))
#         L0 <- L[1, 9:29]
#         x <- 1000
#         L[1, 9:29] <- L0 * x
#         lam <- max(abs(eigen(L, only.values = TRUE)[["values"]]))
#         while (lam <= 1.6) {
#             x <- x * 2
#             L[1, 9:29] <- L0 * x
#             lam <- eigen(L, only.values = TRUE) |>
#                 getElement("values") |> abs() |> max()
#         }
#         return(x)
#     }, .shape, .scale, future.seed = TRUE)
#     return(out)
# }
# max_lam_df <- crossing(shape = seq(0, 5, length.out = 101)[-1],
#                        scale = seq(0, 20, length.out = 101)[-1]) |>
#     mutate(x = get_min_x(shape, scale))
#
# bounds <- c(22, 33, 46, 61, 80)
# # This seems to work okay for categorizing:
# max_lam_df |>
#     mutate(ss = shape * scale,
#            obs_x = case_when(ss > bounds[5] ~ 2^(6-1) * 1000,
#                              ss > bounds[4] ~ 2^(5-1) * 1000,
#                              ss > bounds[3] ~ 2^(4-1) * 1000,
#                              ss > bounds[2] ~ 2^(3-1) * 1000,
#                              ss > bounds[1] ~ 2^(2-1) * 1000,
#                              TRUE ~ 1000)) |>
#     ggplot(aes(obs_x, x)) +
#     geom_abline(slope = 1, intercept = 0, linetype = 2) +
#     geom_point()
# max_lam_df |>
#     mutate(ss = shape * scale,
#            obs_x = case_when(ss > bounds[5] ~ 6,
#                              ss > bounds[4] ~ 5,
#                              ss > bounds[3] ~ 4,
#                              ss > bounds[2] ~ 3,
#                              ss > bounds[1] ~ 2,
#                              TRUE ~ 1000) |>
#            factor()) |>
#     ggplot(aes(ss, x, color = obs_x)) +
#     geom_vline(xintercept = bounds, linetype = 2) +
#     geom_point() +
#     scale_color_viridis_d(begin = 0.2)







get_x <- function(.shape, .scale, .lambda) {
    stopifnot(length(.shape) == length(.scale))
    stopifnot(length(.shape) == length(.lambda))
    stopifnot(all(.shape > 0))
    stopifnot(all(.scale > 0))
    stopifnot(all(.lambda > 0))
    n_sims <- length(.shape)
    p <- progressor(steps = n_sims)
    upper_bounds <- case_when(.shape * .scale > 80 ~ 32000,
                              .shape * .scale > 61 ~ 16000,
                              .shape * .scale > 46 ~ 8000,
                              .shape * .scale > 33 ~ 4000,
                              .shape * .scale > 22 ~ 2000,
                              TRUE ~ 1000)
    out <- future_mapply(\(sh, sc, la, ub) {
        L <- matrix(0, 29, 29)
        L[(row(L) - col(L) == 1)] <- 1
        p_distr_vals <- pweibull(0:21, shape = sh, scale = sc)
        L[1, 9:29] <- diff(p_distr_vals) / sum(diff(p_distr_vals))
        op <- optim_L(L, la, method = "DIRECT_L", max_iters = 1e6L,
                      upper_bound = ub)
        p()
        if (op[3] < 0) return(NA)
        # if (op[2] > 0) return(NA)
        # return(op[1])
        return(tibble(x = op[1], d = op[2]))
    }, .shape, .scale, .lambda, upper_bounds,
    SIMPLIFY = FALSE, future.seed = TRUE)
    return(out)
}


## ****** Uncomment to see plots of x ~ shape * scale  ******
# # Takes ~12 sec with 100 values
# ss_df <- crossing(shape = seq(0, 10, length.out = 10+1)[-1],
#                   scale = seq(0, 20, length.out = 10+1)[-1]) |>
#     mutate(x = get_x(shape, scale, lambda))
#
# ss_df |> filter(is.na(x))
#
#
# ss_df |>
#     mutate(scale = factor(scale)) |>
#     ggplot(aes(log(shape), log(x), color = scale)) +
#     geom_line() +
#     scale_color_viridis_d(begin = 0.2)
#
# ss_df |>
#     mutate(shape = factor(shape)) |>
#     ggplot(aes((scale), (x), color = shape)) +
#     geom_line() +
#     scale_color_viridis_d(begin = 0.2)
#
# ss_df



# ----------------------------------------------------------------------*
# POLYNOMIAL REGRESSION ----
# ----------------------------------------------------------------------*


# This approach works better than `ss_df` above for statistical fitting:
n_ss <- 500L
# Takes ~2.7 min
t0 <- Sys.time()
set.seed(63457)
ss_fit_df <- tibble(shape = runif(n_ss, 0, 5),
                    scale = runif(n_ss, 0, 20),
                    lambda = runif(n_ss, 1.05, 1.6)) |>
    mutate(x = get_x(shape, scale, lambda)) |>
    unnest(x)
t1 <- Sys.time()
t1 - t0; rm(t0, t1)
# should have no rows:
ss_fit_df |> filter(d > 1e-9)

# This is a separate dataset to test the method:
set.seed(1534679)
test_df <- tibble(shape = runif(100, 0, 5),
                  scale = runif(100, 0, 20),
                  lambda = runif(100, 1.05, 1.6)) |>
    mutate(x = get_x(shape, scale, lambda)) |>
    unnest(x)
# should have no rows:
test_df |> filter(d > 1e-9)


# I now test for which polynomial degree to use for each parameter by seeing
# which permutation of each parameter's degree ranging from 1:6
# minimized the sum of the absolute differences between the predicted
# and real `x`.
k_perms <- expand.grid(shape = 1:6, scale = 1:6, lambda = 1:6)
ss <- numeric(nrow(k_perms))
for (k in 1:nrow(k_perms)) {
    # k <- 3
    mod <- lm(log(x) ~ poly(shape, k_perms$shape[[k]], raw = TRUE) *
                  poly(scale, k_perms$scale[[k]], raw = TRUE) *
                  poly(lambda, k_perms$lambda[[k]], raw = TRUE), ss_fit_df)
    x_fit <- suppressWarnings(exp(predict(mod, newdata = test_df[,c("shape", "scale", "lambda")])))
    ss[k] <- sum(abs(x_fit - test_df$x))
}

# This one wins:
k_perms[which(ss == min(ss)),]
#     shape scale lambda
# 101     5     5      3

mod <- lm(log(x) ~ poly(shape, 5, raw = TRUE) *
              poly(scale, 5, raw = TRUE) *
              poly(lambda, 3, raw = TRUE), ss_fit_df)
x_fit <- exp(predict(mod, newdata = test_df[,c("shape", "scale", "lambda")]))

test_df |>
    mutate(x_fit = x_fit) |>
    ggplot(aes(x, x_fit)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "gray70") +
    geom_point()


# Works nicely!



cnames <- names(coef(mod))[!grepl(":|Intercept", names(coef(mod)))] |>
    str_remove_all("\\s+")
if (!all(grepl("^poly", cnames)) || !all(grepl("raw=TRUE", cnames))) {
    stop("all parameters must be polynomials (i.e., use poly(X, k, raw = TRUE)")
}
k_shape <- sum(grepl("shape", cnames))
k_scale <- sum(grepl("scale", cnames))
k_lambda <- sum(grepl("lambda", cnames))

coefs <- coef(mod) |>
    (\(x) {
        max_digits <- ceiling(log10(max(k_shape, k_scale, k_lambda)+1))
        rm_str <- sapply(1:max_digits, \(i) {
            paste0(",", paste(rep(".?", i), collapse = ""), ",raw=TRUE")
        }) |>
            paste(collapse = "|")
        names(x) <- names(x) |>
            str_remove_all("\\s+|\\(|\\)") |>
            str_remove_all(paste0("poly|", rm_str))
        return(x)
    })()


b0 <- coefs[["Intercept"]]

b_shape <- sapply(1:k_shape, \(x) coefs[[paste0("shape",x)]])
b_scale <- sapply(1:k_scale, \(x) coefs[[paste0("scale",x)]])
b_lambda <- sapply(1:k_lambda, \(x) coefs[[paste0("lambda",x)]])

# Two-way interactions:
way2_mats <- list("shape:scale" = matrix(0, k_shape, k_scale),
                  "shape:lambda" = matrix(0, k_shape, k_lambda),
                  "scale:lambda" = matrix(0, k_scale, k_lambda)) |>
    imap(\(imat, pars) {
        pars <- str_split(pars, ":")[[1]]
        for (i in 1:nrow(imat)) for (j in 1:ncol(imat)) {
            imat[i,j] <- coefs[[paste0(pars[1], i, ":", pars[2], j)]]
        }
        return(imat)
    })

# Three-way interactions:
way3_array <- array(0, dim = c(k_shape, k_scale, k_lambda)) |>
    (\(imat) {
        for (i in 1:k_shape) for (j in 1:k_scale) for (k in 1:k_lambda) {
            imat[i,j,k] <- coefs[[sprintf("shape%i:scale%i:lambda%i", i, j, k)]]
        }
        return(imat)
    })()
way3_array[is.na(way3_array)] <- 0


poly_fits <- list("b0" = b0,
                  "b_shape" = b_shape,
                  "b_scale" = b_scale,
                  "b_lambda" = b_lambda,
                  "two_way" = way2_mats,
                  "three_way" = way3_array,
                  "model" = mod)



usethis::use_data(poly_fits, overwrite = TRUE)
