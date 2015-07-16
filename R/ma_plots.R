#' List Count MA
#'
#' Plot MA plots from a list of count matrices
#'
#' @param data.list List of count matrices to plot
#'
#' @return List of ggplot2 object containing MA plots
#'
#' @importFrom magrittr "%>%"
#'
#' @export
listCountMA <- function(data.list) {

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
