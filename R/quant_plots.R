#' List Density
#'
#' Plot densities from a list of matrices
#'
#' @param data.list List of matrices to plot
#'
#' @return List of ggplot2 object containing histograms
#'
#' @importFrom magrittr "%>%"
#'
#' @export
listDensity <- function(data.list) {

    plots <- list()

    for (name in names(data.list)) {

        gg <- countDensity(data.list[[name]])

        switch(
            name,

            "log" = {
                gg <- gg + ggplot2::xlab("log(Counts + 1)")
            },

            "logCPM" = {
                gg <- gg + ggplot2::xlab("log CPM")
            },

            "rlog" = {
                gg <- gg + ggplot2::xlab("rlog Counts")
            },

            "vst" = {
                gg <- gg + ggplot2::xlab("VST Counts")
            }
        )

        plots[[name]] <- gg
    }

    gg <- data.list %>%
          combineMatrices %>%
          magrittr::set_colnames(c("Gene", "Sample", "Counts", "Set")) %>%
          ggplot2::ggplot(aes(x = Counts, fill = Sample, colour = Sample)) +
          ggplot2::geom_density(alpha = 0.3) +
          ggplot2::facet_wrap(~ Set)

    plots[["combined"]] <- gg

    return(plots)
}

#' Count Density
#'
#' Produce a ggplot2 object plotting densities from a matrix of counts
#'
#' @param data Matrix to plot
#'
#' @return ggplot2 object containing the density plot
#'
#' @export
countDensity <- function(data) {

    gg <- data %>%
          lengthenMatrix %>%
          magrittr::set_colnames(c("Gene", "Sample", "Counts")) %>%
          ggplot2::ggplot(aes(x = Counts, fill = Sample, colour = Sample)) +
          ggplot2::geom_density(alpha = 0.3)

    return(gg)
}
