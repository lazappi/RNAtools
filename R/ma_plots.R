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

#' List Results MA
#'
#' Plot MA plots from a list of differential expression results
#'
#' @param data.list List of results to plot
#' @param alpha     Significance level for labelling differentially
#'                  expressed genes
#'
#' @return List of ggplot2 object containing MA plots
#'
#' @importFrom magrittr "%>%"
#'
#' @export
listResultsMA <- function(data.list, alpha = 0.05) {

    plots <- list()

    for (name in names(data.list)) {

        data = data.list[[name]]

        gg <- resultsMA(data, method = name, alpha = alpha)

        plots[[name]] <- gg
    }

    gg <- data.list %>%


    plots[["combined"]] <- gg

    return(plots)
}

#' Results MA
#'
#' Produce an MA plot from results of a differential expression test
#'
#' @param results Results to plot
#' @param method  Method used to produce the results
#' @param alpha   Significance level for labelling differentially
#'                expressed genes
#' @param plot    Boolean, if true return plot, if false return plot data
#'
#' @return ggplot2 object containg MA plot
#'
#' @importFrom magrittr "%>%"
#'
#' @export
resultsMA <- function(results, method = c("edgeR", "DESeq", "DESeq2", "voom"),
                      alpha = 0.05, plot = TRUE) {

    # Check that a valid method has been given
    if (missing(method)) {
        stop("Differential expression method must be specified")
    } else {
        method <- match.arg(method)
    }

    switch(
        method,

        edgeR = {
            plot.data <- results %>%
                         dplyr::select(FoldChange   = logFC,
                                       Abundance    = logCPM,
                                       Significance = FDR) %>%
                         dplyr::mutate(DE = Significance < alpha)

            xlabel <- "logCPM"
            ylabel <- "logFC"
        },

        DESeq = {
            plot.data <- results %>%
                         dplyr::select(FoldChange   = log2FoldChange,
                                       Abundance    = baseMean,
                                       Significance = padj) %>%
                         dplyr::mutate(Abundance = log2(Abundance)) %>%
                         dplyr::mutate(DE = Significance < alpha)

            xlabel <- "log2 mean of normalised counts"
            ylabel <- expression(log[2] ~ fold ~ change)
        },

        DESeq2 = {
            plot.data <- results %>%
                         data.frame %>%
                         dplyr::select(FoldChange   = log2FoldChange,
                                       Abundance    = baseMean,
                                       Significance = padj) %>%
                         dplyr::mutate(Abundance = log2(Abundance)) %>%
                         dplyr::mutate(DE = Significance < alpha)

            xlabel <- "log2 mean of normalised counts"
            ylabel <- expression(log[2] ~ fold ~ change)
        },

        voom = {
            plot.data <- results %>%
                         dplyr::select(FoldChange   = logFC,
                                       Abundance    = AveExpr,
                                       Significance = adj.P.Val) %>%
                         #dplyr::mutate(Abundance = log2(Abundance)) %>%
                         dplyr::mutate(DE = Significance < alpha)

            xlabel <- "log average expression"
            ylabel <- "logFC"
        }
    )

    gg <- ggplot2::ggplot(plot.data,
                          ggplot2::aes(x = Abundance, y = FoldChange,
                                       colour = DE)) +
          ggplot2::geom_rect(aes(xmin = -Inf, xmax = Inf,
                                 ymin = -2,   ymax = 2),
                             fill = "grey", alpha = 0.5, size = 0) +
          ggplot2::geom_point() +
          ggplot2::geom_hline(yintercept = 0, colour = "blue") +
          ggplot2::scale_colour_manual(values = c("black", "red")) +
          #ggplot2::scale_y_continuous(limits = c(-5, 5)) +
          ggplot2::xlab(xlabel) +
          ggplot2::ylab(ylabel) +
          ggplot2::theme(legend.position = "none")

    if (plot) {
        return(gg)
    } else {
        return(plot.data)
    }
}
