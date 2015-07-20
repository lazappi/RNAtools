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
          ggplot2::geom_histogram(binwidth = 0.025) +
          ggplot2::xlab("p-value") +
          ggplot2::ylab("Counts")

    return(gg)
}
