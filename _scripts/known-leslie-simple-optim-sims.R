
#'
#' Simulations to figure out the best combination of parameters to use
#' for `winnowing_optim` function.
#'

source("_scripts/known-leslie-simple-preamble.R")




# =============================================================================*
# functions ----
# =============================================================================*


one_optim_sim <- function(adjust_L, no_K, n_repros, sigma_s, sigma_f, box_optim,
                          n_bevals, n_boxes, seed, overwrite) {

    # This is like match.call but includes default arguments and
    # evaluates everything
    call_ <- mget(names(formals()), sys.frame(sys.nframe()))

    rds_file <- sprintf("_data/optim_sims/%s.rds",
                        paste(adjust_L, no_K, n_repros, sigma_s, sigma_f,
                              box_optim, n_bevals, n_boxes, sep = "_")) |>
        str_replace_all("::", "-")
    # Stop this function if no overwriting desired, file exists, and file
    # reads without error or warning:
    if (!overwrite && file.exists(rds_file)) {
        x <- tryCatch(read_rds(rds_file),
                      error = function(e) NA, warning = function(w) NA)
        if (!isTRUE(is.na(x))) return(invisible(NULL))
    }

    .bof <- eval(parse(text = box_optim))
    .bc <- eval(formals(winnowing_optim)[["box_control"]])
    if (grepl("^nloptr", box_optim)) {
        VERY_SMALL <- max(.Machine$double.eps, .Machine$double.neg.eps)
        b <- formals(fit_known_leslie)[["trans_base"]]
        m <- ifelse(adjust_L == 1L, "max_leslie_x", "max_fecund_x")
        .bc[["lower"]] <- c(rep(VERY_SMALL, 3) + c(1,0,0),
                            -formals(fit_known_leslie)[["max_surv_x"]]) |>
            base::`[`(1:(2L+adjust_L)) |>
            aphidsync:::trans_pars(trans_base = b)
        .bc[["upper"]] <- c(100, 1-VERY_SMALL,
                            formals(fit_known_leslie)[[m]],
                            formals(fit_known_leslie)[["max_surv_x"]]) |>
            base::`[`(1:(2L+adjust_L)) |>
            aphidsync:::trans_pars(trans_base = b)
    }

    set.seed(seed)
    t0 <- Sys.time()
    sim_df <- test_all_sims(300L, LL[[n_repros]],
                            sigma_s = sigma_s,
                            sigma_f = sigma_f,
                            no_K = no_K,
                            adjust_L = adjust_L,
                            optim_args = list(box_optim = .bof,
                                              box_control = .bc,
                                              n_bevals = n_bevals,
                                              n_boxes = n_boxes),
                            show_progress = FALSE)
    t1 <- Sys.time()
    dt <- as.numeric(difftime(t1, t0, units = "min"))

    ## Eventually you'll summarize like this, but since these simulations take
    ## so long, I'm going to save all the data.
    # out <- out |>
    #     group_by(param, sigma_s, sigma_f) |>
    #     summarize(diverg = mean(abs(fit_diffs)),
    #               sd_diverg = sd(abs(fit_diffs)),
    #               accur = mean(abs(obs - fit)),
    #               sd_accur = sd(abs(obs - fit)),
    #               .groups = "drop") |>
    #     mutate(minutes = dt,
    #            adjust_L = adjust_L,
    #            no_K = no_K,
    #            n_repros = n_repros,
    #            box_optim = box_optim,
    #            n_bevals = n_bevals,
    #            n_boxes = n_boxes) |>
    #     select(adjust_L, no_K, n_repros, sigma_s, sigma_f, box_optim,
    #            n_bevals, n_boxes, param, everything())

    out <- list(sims = sim_df, time = dt, call = call_)

    write_rds(out, rds_file)

    invisible(NULL)
}





all_optim_sims <- function(option_df, overwrite) {
    arg_names <- names(formals(one_optim_sim)) |> keep(\(x) x != "overwrite")
    if (!all(arg_names %in% colnames(option_df)) ||
        !all(colnames(option_df) %in% arg_names)) {
        stop(paste("\nERROR: option_df must contain these and only these",
                   "column names:", paste(arg_names, collapse = ", ")))
    }
    stopifnot(is.logical(overwrite) && length(overwrite) == 1)
    option_df[["overwrite"]] <- overwrite
    n_opt_combos <- nrow(option_df)
    pr <- progressor(steps = n_opt_combos)
    for (i in 1:n_opt_combos) {
        do.call(one_optim_sim, as.list(option_df[i,]))
        pr()
    }
    invisible(NULL)
}






# =============================================================================*
# simulations ----
# =============================================================================*


overwrite <- FALSE

# Took ~ 5 hours!
if (length(list.files("_data/optim_sims", "*.rds$")) < 54L) {
    t0 <- Sys.time()
    set.seed(1148372050)
    crossing(adjust_L = 0:2,
             no_K = FALSE,
             n_repros = 2,
             sigma_s = sigma_s_vals[[2]],
             sigma_f = sigma_f_vals[[2]],
             box_optim = c("nloptr::bobyqa", "stats::optim"),
             n_bevals = 100L * c(1L, 2L, 5L),
             n_boxes = 1000L * c(1L, 2L, 5L)) |>
        mutate(seed = sample.int(2^31-1, n())) |>
        all_optim_sims(overwrite = overwrite)
    t1 <- Sys.time()
    cat("==================\n")
    cat("TOTAL TIME\n")
    cat(sprintf("%.2f hours\n", as.numeric(difftime(t1, t0, units = "hours"))))
    cat("==================\n")
    rm(t0, t1)
}




# Takes ~18 sec
opt_sim_df <- list.files("_data/optim_sims", "*.rds$", full.names = TRUE) |>
    map_dfr(\(x) {
        sim_obj <- read_rds(x)
        out <- sim_obj[["sims"]] |>
            split(~ param + sigma_s + sigma_f, drop = TRUE) |>
            map_dfr(\(sx) {
                .diverg <- - log10(abs(sx$fit_diffs))
                .accur <- - log10(abs(sx$obs - sx$fit))
                d_boot <- smean.cl.boot(.diverg, B = 2000)
                a_boot <- smean.cl.boot(.accur, B = 2000)
                sx[1,] |>
                    select(param, sigma_s, sigma_f) |>
                    mutate(diverg = mean(.diverg),
                           lo_diverg = d_boot[["Lower"]],
                           hi_diverg = d_boot[["Upper"]],
                           accur = mean(.accur),
                           lo_accur = a_boot[["Lower"]],
                           hi_accur = a_boot[["Upper"]])
            }) |>
            mutate(minutes = sim_obj[["time"]],
                   adjust_L = sim_obj[["call"]][["adjust_L"]],
                   no_K = sim_obj[["call"]][["no_K"]],
                   n_repros = sim_obj[["call"]][["n_repros"]],
                   box_optim = sim_obj[["call"]][["box_optim"]],
                   n_bevals = sim_obj[["call"]][["n_bevals"]],
                   n_boxes = sim_obj[["call"]][["n_boxes"]]) |>
            select(adjust_L,
                   # no_K, n_repros, sigma_s, sigma_f, ## << these didn't vary in sims
                   box_optim,
                   n_bevals, n_boxes, param, diverg:hi_accur)
        return(out)
    }) |>
    mutate(box_optim = factor(box_optim),
           adjust_L = factor(adjust_L, levels = 0:2,
                             labels = paste0(0:2, "-param scaling")))



.param <- "width99"
.y <- "accur"

opt_sim_df |>
    filter(param == .param) |>
    ggplot(aes(n_boxes, .data[[.y]], color = factor(n_bevals))) +
    geom_point(position = position_dodge(width = 20)) +
    geom_errorbar(aes(ymin = .data[[paste0("lo_", .y)]],
                      ymax = .data[[paste0("hi_", .y)]]),
                  width = 10, position = position_dodge(width = 20)) +
    geom_line() +
    scale_color_viridis_d(begin = 0.2) +
    ggtitle(.param) +
    ylab(ifelse(.y == "diverg", "-- log<sub>10</sub>(Divergence among polished optimizations)",
                "-- log<sub>10</sub>(Absolute difference between fit and observed)")) +
    facet_wrap(~ box_optim + adjust_L) +
    theme(axis.title.y = element_markdown())
