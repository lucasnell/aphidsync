
# Used to check for NA values below:
na_check <- function(x, n) {
    if (length(x) != 1) {
        if (any(is.na(x))) stop(paste0("\n", n, " contains NAs.\n"))
    } else if (length(x) == 1) {
        if (is.na(x)) stop(paste0("\n", n, " is NA.\n"))
    }
}

# Used to check for too low or high values below:
min_max_check <- function(x, n, .min, .max) {
    if (!is.na(.min)) {
        if (any(x < .min)) {
            if (length(x) > 1) {
                stop(paste0("\n", n, " contains values below the minimum ",
                            "allowed (", .min, ").\n"))
            } else {
                stop(paste0("\n", n, " is below the minimum allowed value (",
                            .min, ").\n"))
            }
        }
    }
    if (!is.na(.max)) {
        if (any(x > .max)) {
            if (length(x) > 1) {
                stop(paste0("\n", n, " contains values above the maximum ",
                            "allowed (", .max, ").\n"))
            } else {
                stop(paste0("\n", n, " is above the maximum allowed value (",
                            .max, ").\n"))
            }
        }
    }
}

#' @param x is object
#' @param n is single string indicating name of object (for use in error messages)
#' @param t is single string indicating type of object
#' @param l is single integer indicating length of object or length-n vector
#'     indicating n dimensions. `NA` results in this being ignored.
#'     Defaults to `1L`.
#' @param .min single number indicating minimum allowed value.
#'     `NA` (the default) results in this being ignored.
#' @param .max single number indicating maximum allowed value.
#'     `NA` (the default) results in this being ignored.
#' @noRd
#'
type_checker <- function(x, n, t, l = 1L, .min = NA, .max = NA) {

    # Check for NAs early because it'll influence the special situation below:
    na_check(x, n)

    # Special situation where desired type is integer, but we input a numeric
    # vector or array which often happens and is typically fine.
    # In this case, we want to verify that converting to integer will not
    # change any values, then we convert `t` to whatever class `x` is:
    if (t == "integer" && inherits(x, c("numeric", "array"))) {
        if (any(x != as.integer(x))) {
            stop(paste("\nERROR:", n, "cannot be properly cast as an integer\n"))
        }
        t <- class(x)[[1]]
    }

    if (!inherits(x, t)) stop(paste("\nERROR:", n, "is not a", t, "\n"))
    if (!is.na(l)) {
        if (length(l) == 1 && length(x) != l) {
            stop(paste("\nERROR:", n, "is not of length", l, "\n"))
        }
        if (length(l) > 1 && !identical(dim(x), l)) {
            stop(paste0("\nERROR:", n, " is not of dimensions c(",
                       paste(l, collapse = ", "), ")\n"))
        }
    }
    min_max_check(x, n, .min, .max)
}
