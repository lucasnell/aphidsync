# This loads the ggplot theme used for the figures.

.onLoad <- function(libname, pkgname) {
    # ggplot theme:
    ggplot2::theme_set( ggplot2::`%+replace%`(
        ggplot2::theme_classic(),
        ggplot2::theme(strip.background = ggplot2::element_blank(),
                       strip.text = ggplot2::element_text(size = 10),
                       axis.title = ggplot2::element_text(color = "black", size = 11),
                       axis.text = ggplot2::element_text(color = "black", size = 9),
                       axis.ticks = ggplot2::element_line(color = "black"),
                       legend.background = ggplot2::element_blank(),
                       plot.title = ggplot2::element_text(size = 14, hjust = 0.5))))
}
