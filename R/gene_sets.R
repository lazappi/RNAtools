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
          ggplot2::coord_flip() +
          ggplot2::guides(shape  = ggplot2::guide_legend(title = "Method"),
                          colour = ggplot2::guide_colourbar(title = "-log Sig")) +
          ggplot2::scale_color_gradient(low = "#3f007d", high = "#ef3b2c") +
          ggplot2::ylab("log Fold Change")

    return(gg)

}

#' Gene Set Table
#'
#' Produce a table showing how many genes in a list of sets are differentially
#' expressed, up-regulated and down-regulated.
#'
#' @param data.list List of differential expression results
#' @param de.set    Set of differentially expressed genes
#' @param gene.sets Named list of gene sets to include in table
#'
#' @return data.frame containing the summary table
#'
#' @importFrom magrittr "%>%"
#'
#' @export
setTable <- function(data.list, de.set, gene.sets) {

    gene.summary <- geneSummary(data.list)

    genes.up <-  gene.summary %>%
                 dplyr::filter(meanFC > 0) %>%
                 as.data.frame() %>%
                 magrittr::extract(, "Gene")

    genes.dn <-  gene.summary %>%
                 dplyr::filter(meanFC <= 0) %>%
                 as.data.frame() %>%
                 magrittr::extract(, "Gene")

    de.set.up <- gene.summary %>%
                 dplyr::filter(Gene %in% de.set) %>%
                 dplyr::filter(meanFC > 0) %>%
                 as.data.frame() %>%
                 magrittr::extract(, "Gene")

    de.set.dn <- gene.summary %>%
                 dplyr::filter(Gene %in% de.set) %>%
                 dplyr::filter(meanFC <= 0) %>%
                 as.data.frame() %>%
                 magrittr::extract(, "Gene")

    table.data <- lapply(gene.sets,
                         function(set) {
                             total <- length(set)
                             up    <- length(intersect(set, genes.up))
                             dn    <- length(intersect(set, genes.dn))
                             de    <- length(intersect(set, de.set))
                             de.up <- length(intersect(set, de.set.up))
                             de.dn <- length(intersect(set, de.set.dn))
                             return(c(total, up, dn, de, de.up, de.dn))
                         }
                    )

    table.data <- table.data %>% unlist %>% matrix(ncol = 6, byrow = TRUE)

    table <- data.frame(Set      = names(gene.sets),
                        NumGenes = table.data[, 1],
                        Up       = table.data[, 2],
                        Down     = table.data[, 3],
                        DE       = table.data[, 4],
                        DEUp     = table.data[, 5],
                        DEDown   = table.data[, 6])
    return(table)
}
