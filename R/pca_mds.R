#' Plot PCA
#'
#' Produce a PCA plot from a matrix using ggplot2
#'
#' @param data   Matrix of data to plot
#' @param top    Number of rows with highes variance to select for plotting
#' @param groups Vector of groups assigned to sample columns
#'
#' @return ggplot2 object containg the PCA plot
#'
#' @importFrom magrittr "%>%"
#'
#' @export
plotPCA <- function(data, top = nrow(data), groups = colnames(data)) {

    top.data <- data %>%
                data.frame %>%
                dplyr::mutate(var = genefilter::rowVars(data)) %>%
                dplyr::top_n(top, var) %>%
                dplyr::select(-var)

    PCA.data <- prcomp(t(top.data), scale = FALSE)

    percent.var <- round(100 * PCA.data$sdev ^ 2 / sum(PCA.data$sdev ^ 2), 1)

    print(PCA.data$sdev)

    plot.data <- data.frame(PC1    = PCA.data$x[, 1],
                            PC2    = PCA.data$x[, 2],
                            Sample = colnames(data),
                            Group  = groups)

    gg <- ggplot2::ggplot(plot.data,
                          ggplot2::aes(x = PC1, y = PC2, colour = Group,
                                       label = Sample)) +
          ggplot2::geom_text() +
          ggplot2::ggtitle(paste("PC1 vs PC2, top", top, "variable genes")) +
          ggplot2::xlab(paste0("PC1, VarExp: ", round(percent.var[1], 4), "%")) +
          ggplot2::ylab(paste0("PC2, VarExp: ", round(percent.var[2], 4), "%")) +
          ggplot2::theme(axis.title = ggplot2::element_text(size = 20),
                         axis.text  = ggplot2::element_text(size = 15),
                         plot.title = ggplot2::element_text(size = 30,
                                                            face = "bold"),
                         legend.text = ggplot2::element_text(size = 15),
                         legend.title = ggplot2::element_text(size = 15))

    return(gg)
}
