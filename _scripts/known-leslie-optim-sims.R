
#'
#' Simulations to figure out the best combination of parameters to use
#' for `winnowing_optim` function, for more complex Leslie matrices.
#'

source("_scripts/known-leslie-simple-preamble.R")

sim_rds_dir <- "_data/optim_sims_complex"
opt_sim_summ_rds <- "_data/opt_sim_df_complex.rds"




# =============================================================================*
# functions ----
# =============================================================================*


one_optim_sim <- function(adjust_L, box_optim, n_bevals,
                          n_stages, n_repros,
                          seed, overwrite) {

    # adjust_L = 1L; box_optim <- fine_optim <- polished_optim <- "stats::optim";
    # n_bevals = 100; seed = 46735890; overwrite = FALSE
    # rm(adjust_L, box_optim, fine_optim, polished_optim, n_bevals, seed,
    #    overwrite, call_, rds_file, mk_cntrl, conts, funs, t0, t1, sim_df,
    #    dt, out)

    # This is like match.call but includes default arguments and
    # evaluates everything
    call_ <- mget(names(formals()), sys.frame(sys.nframe()))

    stopifnot(n_stages >= n_repros)

    rds_file <- sprintf("_data/optim_sims_complex/%s.rds",
                        paste(adjust_L, box_optim, n_bevals,
                              n_stages, n_repros, sep = "-")) |>
        str_replace_all("::", "_")
    # Stop this function if no overwriting desired, file exists, and file
    # reads without error or warning:
    if (!overwrite && file.exists(rds_file)) {
        x <- tryCatch(read_rds(rds_file),
                      error = function(e) NA, warning = function(w) NA)
        if (!isTRUE(is.na(x))) return(invisible(NULL))
    }


    # Make control list for box optimization:
    box_cntrl <- eval(formals(winnowing_optim)[["box_control"]])
    if (grepl("^nloptr|^minqa", box_optim)) {
        VERY_SMALL <- max(.Machine$double.eps, .Machine$double.neg.eps)
        tb <- formals(fit_known_leslie)[["trans_base"]]
        msx <- formals(fit_known_leslie)[["max_surv_x"]]
        mfx <- formals(fit_known_leslie)[["max_fecund_x"]]
        box_cntrl[["lower"]] <- c(rep(VERY_SMALL, 3) + c(1,0,0), -msx) |>
            base::`[`(1:(2L+adjust_L)) |>
            aphidsync:::trans_pars(trans_base = tb)
        box_cntrl[["upper"]] <- c(100, 1-VERY_SMALL, mfx, msx) |>
            base::`[`(1:(2L+adjust_L)) |>
            aphidsync:::trans_pars(trans_base = tb)
    }

    # Optimizer function:
    box_fun <- eval(parse(text = box_optim))

    set.seed(seed)
    t0 <- Sys.time()
    sim_df <- test_all_sims(300L, L = NULL,
                            n_stages = n_stages,
                            n_repros = n_repros,
                            sigma_s = sigma_s_vals[[2]],
                            sigma_f = sigma_f_vals[[2]],
                            no_K = FALSE,
                            adjust_L = adjust_L,
                            optim_args = list(box_optim = box_fun,
                                              box_control = box_cntrl,
                                              n_bevals = n_bevals),
                            show_progress = FALSE)
    t1 <- Sys.time()
    dt <- as.numeric(difftime(t1, t0, units = "min"))

    out <- list(sims = sim_df, minutes = dt, call = call_)

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

# Takes ~ hours!
if (overwrite || length(list.files(sim_rds_dir, "*.rds$")) < 198L) {
    t0 <- Sys.time()
    set.seed(1607759038)
    crossing(adjust_L = 0:2,
             box_optim = c("minqa::bobyqa", "stats::optim"),
             n_bevals = 100L * c(1L, 5L, 10L),
             n_stages = c(5, 10, 20),
             n_repros = c(1, 2, 5, 10)) |>
        filter(n_stages >= n_repros) |>
        mutate(seed = sample.int(2^31-1, n())) |>
        all_optim_sims(overwrite = overwrite)
    t1 <- Sys.time()
    cat("==================\n")
    cat("TOTAL TIME\n")
    cat(sprintf("%.2f hours\n", as.numeric(difftime(t1, t0, units = "hours"))))
    cat("==================\n")
    rm(t0, t1)
}


if (overwrite || !file.exists(opt_sim_summ_rds)) {
    # Takes ~18 sec
    opt_sim_df <- list.files(sim_rds_dir, "*.rds$", full.names = TRUE) |>
        map_dfr(\(x) {

            sim_obj <- read_rds(x)

            arg_names <- names(sim_obj[["call"]]) |>
                discard(\(x) x %in% c("seed", "overwrite"))

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
                })

            for (an in arg_names) {
                out[[an]] <- sim_obj[["call"]][[an]]
            }

            out <- out |>
                select(all_of(arg_names), param, diverg:hi_accur)

            time_row <- bind_cols(out[1,arg_names],
                                  tibble(param = "minutes",
                                         diverg = sim_obj[["minutes"]],
                                         accur = sim_obj[["minutes"]]))
            out <- bind_rows(out, time_row)
            return(out)
        }) |>
        mutate(across(ends_with("_optim"),
                      \(x) factor(x, levels = c("minqa::bobyqa", "stats::optim"),
                                  labels = c("bobyqa", "optim"))),
               adjust_L = factor(adjust_L, levels = 0:2,
                                 labels = paste0(0:2, "-param")))
    write_rds(opt_sim_df, opt_sim_summ_rds)
} else {
    opt_sim_df <- read_rds(opt_sim_summ_rds)
}





optim_plotter <- function(.p, .y, .facet_nrow = 1L, .y_lims = NULL) {

    # `.y` doesn't matter in this case:
    if (missing(.y) && .p == "minutes") .y <- "diverg"

    param_map <- list(shape = "Initial Beta shape",
                      offset = "Initial Beta offset",
                      K = "Density dependence",
                      width99 = "Width of 99% quantile",
                      med_age = "Median starting age",
                      minutes = "Minutes elapsed")
    ylab_fun <- function(.y, .p) {
        if (.p == "minutes") return("Total elapsed minutes")
        ym <- sprintf("-- log<sub>10</sub>(%s)",
                      c("Divergence among polished optimizations",
                        "Absolute difference between fit and observed")) |>
            set_names(c("diverg", "accur")) |>
            as.list()
        return(ym[[.y]])
    }

    stopifnot(.p %in% names(param_map))
    stopifnot(.y %in% c("diverg", "accur"))

    opt_sim_df |>
        filter(param == .p) |>
        mutate(n_bevals = factor(n_bevals),
               id = interaction(box_optim, fine_optim, polished_optim,
                                drop = TRUE, sep = "\n> ")) |>
        ggplot(aes(n_bevals, .data[[.y]], color = adjust_L)) +
        geom_hline(aes(yintercept = mean(.data[[.y]])),
                   linetype = "22", color = "gray70") +
        geom_errorbar(aes(ymin = .data[[paste0("lo_", .y)]],
                          ymax = .data[[paste0("hi_", .y)]]),
                      width = 0.2, position = position_dodge(width = 0.5),
                      na.rm = TRUE) +
        geom_line(aes(x = as.numeric(n_bevals)),
                  position = position_dodge(width = 0.5)) +
        geom_point(position = position_dodge(width = 0.5)) +
        scale_color_viridis_d("Leslie\nscaling", begin = 0.2) +
        ggtitle(param_map[[.p]]) +
        xlab("Number of evaluations per box") +
        scale_y_continuous(ylab_fun(.y, .p), limits = .y_lims) +
        facet_wrap(~ id, nrow = .facet_nrow, scales = "fixed") +
        theme(axis.title.y = element_markdown(),
              plot.title = element_text(size = 16, face = "bold"),
              strip.text = element_text(size = 11, lineheight = unit(0.8, "lines")))
}


optim_plotter(.p = "minutes", .y_lims = c(0, NA))

optim_plotter(.p = "width99", .y = "diverg")
optim_plotter(.p = "offset", .y = "diverg")

optim_plotter(.p = "width99", .y = "accur")
optim_plotter(.p = "offset", .y = "accur")


# From this, I think that the 2-parameter scaling with bobyqa > optim > optim


.p = "offset"
.y = "accur"

opt_sim_df |>
    filter(param == .p) |>
    mutate(n_bevals = factor(n_bevals),
           id = interaction(box_optim, fine_optim, polished_optim,
                            drop = TRUE, sep = " > ")) |>
    filter(polished_optim == "optim", fine_optim == "optim",
           adjust_L == "2-param") |>
    ggplot(aes(n_bevals, .data[[.y]], color = id)) +
    geom_errorbar(aes(ymin = .data[[paste0("lo_", .y)]],
                      ymax = .data[[paste0("hi_", .y)]]),
                  width = 0.2, na.rm = TRUE) +
    geom_line(aes(x = as.numeric(n_bevals))) +
    geom_point() +
    scale_color_viridis_d(NULL, begin = 0.2) +
    ggtitle(.p) +
    xlab("Number of evaluations per box") +
    scale_y_continuous(.y) +
    theme(axis.title.y = element_markdown(),
          plot.title = element_text(size = 16, face = "bold"),
          legend.position = "top")
