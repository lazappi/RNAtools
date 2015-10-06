#' List p-value plots
#'
#' Plot histograms of p-values from a list of regularised differential
#' expression results
#'
#' @param data.list List of results to plot
#' @param alpha     Significance level for shading
#'
#' @return List of ggplot2 objects containing p-value histograms
#'
#' @importFrom magrittr "%>%"
#'
#' @export
listPlotPvals <- function(data.list, alpha = 0.05) {

    plots <- list()

    # Produce individual plots
    for (name in names(data.list)) {

        data <- data.list[[name]]

        gg <- plotPvals(data, method = "regular")

        plots[[name]] <- gg
    }

    # Produce combine plot
    gg <- data.list %>%
          combineMatrices(lengthen = FALSE) %>%
          ggplot2::ggplot(ggplot2::aes(x = pValue)) +
          ggplot2::annotate("rect",
                            xmin = -Inf, xmax = alpha,
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
#' @param alpha   Significance level for shading
#'
#' @return ggplot2 object containg p-value histogram
#'
#' @importFrom magrittr "%>%"
#'
#' @export
plotPvals <- function(results,
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

    gg <- ggplot2::ggplot(results, ggplot2::aes(x = pValue)) +
          ggplot2::annotate("rect",
                            xmin = -Inf, xmax = alpha,
                            ymin = -Inf, ymax = Inf,
                            alpha = 0.2, fill = "red", size = 0) +
          ggplot2::geom_histogram(binwidth = 0.025) +
          ggplot2::xlab("p-value") +
          ggplot2::ylab("Counts")

    return(gg)
}
