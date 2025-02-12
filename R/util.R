
#' Width of 99th quartile (offset ignored)
#' @noRd
#' @export
width99 <- function(shape, offset) {
    sapply(shape, \(sh) {
        q99 <- qbeta(c(0.005, 0.995), sh, sh)
        return(abs(diff(q99)))
    })
}
#' Proportion of population that's reproducing
#' @noRd
#' @export
p_repro <- function(shape, offset) {
    mapply(\(sh, of) {
        a0 <- beta_starts(sh, of, 1)
        return(sum(a0[10:29]))
    }, sh = shape, of = offset)
}
#' Median age
#' @noRd
#' @export
med_age <- function(shape, offset) {
    mapply(\(sh, of) {
        cs_a0 <- cumsum(beta_starts(sh, of, 1))
        if (any(cs_a0 == 0.5)) {
            idx <- mean(which(cs_a0 == 0.5) + 0:1)
        } else idx <- which(cs_a0 > 0.5)[[1]]
        return(idx)
    }, sh = shape, of = offset)
}
