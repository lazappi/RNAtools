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
          ggplot2::ggplot(aes(x = Counts, fill = Sample, colour = Sample)) +
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
          ggplot2::ylab("")

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
        ggplot2::xlab("")

    return(gg)
}

#' List MA
#'
#' Plot MA plots from a list of matrices
#'
#' @param data.list List of matrices to plot
#'
#' @return List of ggplot2 object containing MA plots
#'
#' @importFrom magrittr "%>%"
#'
#' @export
listMA <- function(data.list) {

    plots <- list()

    for (name in names(data.list)) {

        gg <- countMA(data.list[[name]])

        switch(
            name,

            "log" = {
                gg <- gg + ggplot2::ylab(
                                paste("Average",
                                      expression(log[2](Counts + 1))))
            },

            "logCPM" = {
                gg <- gg + ggplot2::ylab("Average log CPM")
            },

            "rlog" = {
                gg <- gg + ggplot2::ylab("Average rlog Counts")
            },

            "vst" = {
                gg <- gg + ggplot2::ylab("Average VST Counts")
            }
        )

        plots[[name]] <- gg
    }

    gg <- data.list %>%
          lapply(getMAData) %>%
          combineMatrices(lengthen = FALSE) %>%
          magrittr::set_colnames(c("Gene", "Pair", "M", "A", "Set")) %>%
          ggplot2::ggplot(aes(x = A, y = M)) +
          ggplot2::geom_point() +
          ggplot2::geom_hline(color = "blue3") +
          ggplot2::stat_smooth(se = FALSE, method = "loess", color = "red3") +
          ggplot2::facet_wrap(~ Pair + Set) +
          ggplot2::ylab("Difference") +
          ggplot2::xlab("Average Counts")

    plots[["combined"]] <- gg

    return(plots)
}

#' Count MA
#'
#' Produce a ggplot2 object containing an MA plot from a matrix of counts
#'
#' @param data Matrix to plot
#'
#' @return ggplot2 object containing the MA plot
#'
#' @importFrom magrittr "%>%"
#'
#' @export
countMA <- function(data) {

    plot.data <- getMAData(data)

    gg <- plot.data %>%
          ggplot2::ggplot(aes(x = A, y = M)) +
          ggplot2::geom_point() +
          ggplot2::geom_hline(color = "blue3") +
          ggplot2::stat_smooth(se = FALSE, method = "loess", color = "red3") +
          ggplot2::facet_wrap(~ Pair) +
          ggplot2::ylab("Difference") +
          ggplot2::xlab("Average Counts")

    return(gg)
}


getMAData <- function(data) {

    MA.idx <- t(combn(1:ncol(data), 2))

    ma.data <- data.frame(Gene = character(),
                          Pair = character(),
                          M    = numeric(),
                          A    = numeric())

    for (i in seq_along(MA.idx[, 1])) {

        idx1 <- MA.idx[i, 1]
        idx2 <- MA.idx[i, 2]

        name1 <- colnames(data)[idx1]
        name2 <- colnames(data)[idx2]

        M <-  data[, idx1] - data[, idx2]
        A <- (data[, idx1] + data[, idx2]) / 2

        Pair <- rep(paste(name1, "vs", name2), nrow(data))

        ma.data <- rbind(ma.data,
                         data.frame(Gene = rownames(data), Pair, M, A))

    }

    return(ma.data)
}
