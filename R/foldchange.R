#' Plot Fold Change
#'
#' Produces a scatter plot of estimated fold change produced by multiple
#' methods for a selected set of genes. Point shapes indicate method and colour
#' shows significance (red indicating low pvalues) with stars indicating mean
#' fold change.
#'
#' @param data.list List of results to plot
#' @param gene.set  Vector of genes to select for plotting
#'
#' @return ggplot2 object containing the scatter plot
#'
#' @importFrom magrittr "%>%"
#'
#' @export
plotFoldChange <- function(data.list, gene.set) {

    gene.set <- sort(gene.set)

    regular.data.list <- list()

    for (name in names(data.list)) {

        regular.data <- data.list[[name]] %>%
                        regulariseResults(name) %>%
                        dplyr::filter(Gene %in% gene.set)

        regular.data.list[[name]] <- regular.data
    }

    plot.data <- combineMatrices(regular.data.list, lengthen = FALSE)

    order <- plot.data %>%
             dplyr::group_by(Gene) %>%
             dplyr::summarise(meanFC = mean(FoldChange)) %>%
             as.data.frame() %>%
             magrittr::extract(, "meanFC") %>%
             order()

    addMean <- function() {
        ggplot2::stat_summary(fun.y = mean, colour = "orange", geom = "point",
                              size = 12, shape = "*")
    }

    gg <- plot.data %>%
          dplyr::mutate(Gene = factor(Gene, levels = gene.set[order])) %>%
          ggplot2::ggplot(ggplot2::aes(x = Gene,
                                       y = FoldChange,
                                       colour = -log(Significance),
                                       shape = matrix)) +
          ggplot2::geom_point(size = 4) +
          ggplot2::geom_hline(ggplot2::aes(yintercept = 0), colour = "blue") +
          addMean() +
          ggplot2::annotate("rect",
                            ymin = -2,     ymax = 2,
                            xmin = -Inf,   xmax = Inf,
                            fill = "blue", alpha = 0.1, size = 0) +
          ggplot2::theme_bw() +
          ggplot2::coord_flip() +
          ggplot2::guides(shape = ggplot2::guide_legend(title = "Method"),
                          colour = "none") +
          ggplot2::scale_color_gradient(low = "#3f007d", high = "#ef3b2c") +
          ggplot2:::theme(
              axis.title   = ggplot2::element_text(size   = 20,
                                                   face   = "bold"),
              axis.text    = ggplot2::element_text(size   = 12,
                                                   colour = "grey50"),
              axis.text.y  = ggplot2::element_text(size   = 8),
              plot.title   = ggplot2::element_text(size   = 40,
                                                   face   = "bold",
                                                   vjust  = 1),
              panel.grid   = ggplot2::element_blank(),
              panel.border = ggplot2::element_blank(),
              axis.ticks   = ggplot2::element_blank(),
              legend.position = c(0.85, 0.2),
              legend.key.size = grid::unit(1, "cm"),
              legend.title    = ggplot2::element_text(size = 12),
              legend.text     = ggplot2::element_text(size = 12))

    return(gg)

}
