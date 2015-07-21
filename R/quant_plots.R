#' List Density
#'
#' Plot densities from a list of matrices
#'
#' @param data.list List of matrices to plot
#'
#' @return List of ggplot2 object containing density plots
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
                gg <- gg + ggplot2::xlab(expression(log[2](Counts + 1)))
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
#' @importFrom magrittr "%>%"
#'
#' @export
countDensity <- function(data) {

    gg <- data %>%
          lengthenMatrix %>%
          magrittr::set_colnames(c("Gene", "Sample", "Counts")) %>%
          ggplot2::ggplot(ggplot2::aes(x = Counts, fill = Sample,
                                       colour = Sample)) +
          ggplot2::geom_density(alpha = 0.3)

    return(gg)
}

#' List Boxplots
#'
#' Plot boxplots from a list of matrices
#'
#' @param data.list List of matrices to plot
#'
#' @return List of ggplot2 object containing boxplots
#'
#' @importFrom magrittr "%>%"
#'
#' @export
listBoxplots <- function(data.list) {

    plots <- list()

    for (name in names(data.list)) {

        gg <- countBoxplots(data.list[[name]])

        switch(
            name,

            "log" = {
                gg <- gg + ggplot2::ylab(expression(log[2](Counts + 1)))
            },

            "logCPM" = {
                gg <- gg + ggplot2::ylab("log CPM")
            },

            "rlog" = {
                gg <- gg + ggplot2::ylab("rlog Counts")
            },

            "vst" = {
                gg <- gg + ggplot2::ylab("VST Counts")
            }
        )

        plots[[name]] <- gg
    }

    gg <- data.list %>%
          combineMatrices %>%
          magrittr::set_colnames(c("Gene", "Sample", "Counts", "Set")) %>%
          ggplot2::ggplot(aes(x = Sample, y = Counts, fill = Sample)) +
          ggplot2::geom_boxplot() +
          ggplot2::facet_wrap(~ Set) +
          ggplot2::ylab("") +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                             hjust =  1))

    plots[["combined"]] <- gg

    return(plots)
}

#' Count Boxplot
#'
#' Produce a ggplot2 object plotting boxplots from a matrix of counts
#'
#' @param data Matrix to plot
#'
#' @return ggplot2 object containing the boxplots
#'
#' @importFrom magrittr "%>%"
#'
#' @export
countBoxplots <- function(data) {

    gg <- data %>%
        lengthenMatrix %>%
        magrittr::set_colnames(c("Gene", "Sample", "Counts")) %>%
        ggplot2::ggplot(aes(x = Sample, y = Counts, fill = Sample)) +
        ggplot2::geom_boxplot() +
        ggplot2::xlab("") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                           hjust =  1))

    return(gg)
}


