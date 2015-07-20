#' List p-value plots
#'
#' Plot hitograms of p-values from a list of differential expression results
#'
#' @param data.list List of results to plot
#'
#' @return List of ggplot2 objects containing p-value histograms
#'
#' @importFrom magrittr "%>%"
#'
#' @export
listPlotPvals <- function(data.list) {

    plots <- list()
    regular.data.list <- list()

    for (name in names(data.list)) {

        data <- data.list[[name]]

        regular.data <- regulariseResults(data, name)

        gg <- plotPvals(regular.data, method = "regular")

        regular.data.list[[name]] <- regular.data
        plots[[name]] <- gg
    }

    gg <- regular.data.list %>%
          combineMatrices(lengthen = FALSE) %>%
          ggplot2::ggplot(ggplot2::aes(x = pValue)) +
          ggplot2::annotate("rect",
                            xmin = -Inf, xmax = 0.05,
                            ymin = -Inf, ymax = Inf,
                            alpha = 0.2, fill = "red", size = 0) +
          ggplot2::geom_histogram(binwidth = 0.025) +
          ggplot2::facet_wrap(~ matrix) +
          ggplot2::xlab("p-value") +
          ggplot2::ylab("Counts")

    plots[["combined"]] <- gg

    return(plots)
}

#' Plot p-values
#'
#' Produce a histogram of p-values from differential expression results
#'
#' @param results Results to plot
#' @param method  Method used to produce the results
#'
#' @return ggplot2 object containg p-value histogram
#'
#' @importFrom magrittr "%>%"
#'
#' @export
plotPvals <- function(results,
                        method = c("edgeR", "DESeq", "DESeq2",
                                   "voom", "regular")) {

    # Check that a valid method has been given
    if (missing(method)) {
        stop("Method used to produce results must be specified")
    } else {
        method <- match.arg(method)
    }

    if (method != "regular") {
        results <- results %>% regulariseResults(method = method)
    }

    gg <- ggplot2::ggplot(results, ggplot2::aes(x = pValue)) +
          ggplot2::annotate("rect",
                            xmin = -Inf, xmax = 0.05,
                            ymin = -Inf, ymax = Inf,
                            alpha = 0.2, fill = "red", size = 0) +
          ggplot2::geom_histogram(binwidth = 0.025) +
          ggplot2::xlab("p-value") +
          ggplot2::ylab("Counts")

    return(gg)
}
