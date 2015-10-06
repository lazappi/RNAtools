#' List Volcano
#'
#' Plot volcano plots from a list of differential expression results
#'
#' @param data.list List of results to plot
#' @param alpha     Significance level for shading
#'
#' @return List of ggplot2 object containing volcano plots
#'
#' @importFrom magrittr "%>%"
#'
#' @export
listVolcano <- function(data.list, alpha = 0.05) {

    plots <- list()
    regular.data.list <- list()

    # Produce individual plots
    for (name in names(data.list)) {

        data <- data.list[[name]]

        regular.data <- regulariseResults(data, name)

        gg <- volcanoPlot(regular.data, method = "regular")

        regular.data.list[[name]] <- regular.data
        plots[[name]] <- gg
    }

    # Produce combined plot
    gg <- regular.data.list %>%
          combineMatrices(lengthen = FALSE) %>%
          dplyr::mutate(logSig = -log(Significance)) %>%
          dplyr::mutate(logSig = replace(logSig, is.na(logSig), -Inf)) %>%
          dplyr::mutate(absFC = abs(FoldChange)) %>%
          dplyr::mutate(absFC = replace(absFC, !is.finite(absFC),
                                        max(absFC[is.finite(absFC)]))) %>%
          ggplot2::ggplot(ggplot2::aes(x      = FoldChange, y = logSig,
                                       colour = logSig, alpha = absFC)) +
          ggplot2::annotate("rect",
                            xmin = -Inf, xmax = Inf,
                            ymin = -Inf, ymax = -log(alpha),
                            alpha = 0.2, fill = "red", size = 0) +
          ggplot2::annotate("rect",
                            xmin = -2,   xmax = 2,
                            ymin = -Inf, ymax = Inf,
                            alpha = 0.2, fill = "blue", size = 0) +
          ggplot2::geom_point(size = 4) +
          ggplot2::facet_wrap(~ matrix) +
          ggplot2::scale_color_gradient(low = "#3f007d", high = "#ef3b2c") +
          ggplot2::xlab("log Fold Change") +
          ggplot2::ylab("-log Significance") +
          ggplot2::theme(legend.position = "none")


    plots[["combined"]] <- gg

    return(plots)
}

#' Volcano Plot
#'
#' Produce a volcano plot from differential expression results
#'
#' @param results Results to plot
#' @param method  Method used to produce the results
#' @param alpha   Significance level for shading
#'
#' @return ggplot2 object containg volcano plot
#'
#' @importFrom magrittr "%>%"
#'
#' @export
volcanoPlot <- function(results,
                        method = c("edgeR", "DESeq", "DESeq2",
                                   "voom", "regular"),
                        alpha = 0.05) {

    # Check that a valid method has been given
    if (missing(method)) {
        stop("Method used to produce results must be specified")
    } else {
        method <- match.arg(method)
    }

    if (method != "regular") {
        results <- results %>% regulariseResults(method = method)
    }

    # Reshape data from plotting
    plot.data <- results %>%
                 dplyr::mutate(logSig = -log(Significance)) %>%
                 dplyr::mutate(logSig = replace(logSig, is.na(logSig), -Inf)) %>%
                 dplyr::mutate(absFC = abs(FoldChange)) %>%
                 dplyr::mutate(absFC = replace(absFC, !is.finite(absFC),
                                               max(absFC[is.finite(absFC)])))

    gg <- ggplot2::ggplot(plot.data,
                          ggplot2::aes(x      = FoldChange, y = logSig,
                                       colour = logSig, alpha = absFC)) +
          ggplot2::annotate("rect",
                            xmin = -Inf, xmax = Inf,
                            ymin = -Inf, ymax = -log(alpha),
                            alpha = 0.2, fill = "red", size = 0) +
          ggplot2::annotate("rect",
                            xmin = -2,   xmax = 2,
                            ymin = -Inf, ymax = Inf,
                            alpha = 0.2, fill = "blue", size = 0) +
          ggplot2::geom_point(size = 4) +
          ggplot2::scale_color_gradient(low = "#3f007d", high = "#ef3b2c") +
          ggplot2::xlab("log Fold Change") +
          ggplot2::ylab("-log Significance") +
          ggplot2::theme(legend.position = "none")

    return(gg)
}
