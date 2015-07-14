#' List Heatmaps
#'
#' Consturct heatmaps from a list of matrices
#'
#' @param data.list List of matrices to plot
#' @param groups Vector of groups for samples
#'
#' @return List of lists containing objects for heatmaps
#'
#' @export
listHeatmaps <- function(data.list, groups = colnames(data.list[[1]])) {

    plots <- list()

    for (name in names(data.list)) {

        heatmap <- countHeatmap(data.list[[name]], groups = groups)

        plots[[name]] <- heatmap
    }

    return(plots)
}

#' Count heatmap
#'
#' Produce a clustered heatmap from a count matrix
#'
#' @param data Matrix to produce the heatmap from
#' @param groups Vector of groups for samples
#'
#' @return List counting elements of heatmap
#'
#' @export
countHeatmap <- function(data, groups = colnames(data)) {

    dists <- as.matrix(dist(t(data)))
    colnames(dists) <- groups

    heatmap <- makeHeatmap(dists, dist.mat = TRUE)

    return(heatmap)

}

#' Heatmap Dendrogram
#'
#' Produce a dendrogram for a heatmap using ggplot2 and ggdendro
#'
#' @param data Dendrogram data to plot
#' @param rows Boolean, whether the plot is for rows
#'
#' @return ggplot2 object containing the dendrogram
heatmapDendro <- function(data, rows = TRUE) {

    if (rows) {
        gg <- ggdendro::ggdendrogram(data, labels = FALSE, rotate = TRUE) +
              ggplot2::theme(axis.text.x = ggplot2::element_blank())
    } else {
        gg <- ggdendro::ggdendrogram(data, labels = FALSE, rotate = FALSE) +
              ggplot2::theme(axis.text.y = ggplot2::element_blank())
    }

    return(gg)
}

#' Make Heatmap
#'
#' Produce the group of plots required for a heatmap with cluster dendrograms
#' from a matrix
#'
#' @param data     Matrix of data to plot
#' @param dist.mat Boolean, whether the matrix already contains distances
#'
#' @return list of components of the cluster heatmap
#'
#' @importFrom magrittr "%>%"
#'
#' @export
makeHeatmap <- function(data, dist.mat = FALSE) {

    colours <- rev(colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(255))

    if (dist.mat) {
        row.hc <- hclust(as.dist(data),    "ward.D")
        col.hc <- hclust(as.dist(t(data)), "ward.D")
    } else {
        row.hc <- hclust(dist(data),    "ward.D")
        col.hc <- hclust(dist(t(data)), "ward.D")
    }

    row.dendro <- ggdendro::dendro_data.dendrogram(as.dendrogram(row.hc),
                                                   type = "rectangle")
    col.dendro <- ggdendro::dendro_data.dendrogram(as.dendrogram(col.hc),
                                                   type = "rectangle")

    row.plot <- heatmapDendro(row.dendro, rows = TRUE) +
                ggplot2::theme(plot.margin = grid::unit(c(0, 0, 0, 0), "lines"))
    col.plot <- heatmapDendro(col.dendro, rows = FALSE) +
                ggplot2::theme(plot.margin = grid::unit(c(0, 0, 0, 0), "lines"))

    row.ord <- match(row.dendro$labels$label, rownames(data))
    if (dist.mat) {
        col.ord <- row.ord
    } else {
        col.ord <- match(col.dendro$labels$label, colnames(data))
    }

    plot.data <- data %>%
                 data.frame %>%
                 dplyr::select(col.ord) %>%
                 dplyr::slice(row.ord) %>%
                 magrittr::set_rownames(rownames(data)[row.ord]) %>%
                 dplyr::mutate(rowname =
                                   factor(rownames(data)[row.ord],
                                          levels = rownames(data)[row.ord])) %>%
                 tidyr::gather(key = rowname) %>%
                 magrittr::set_colnames(c("row", "col", "value"))

    gg <- ggplot2::ggplot(plot.data, ggplot2::aes(row, col)) +
          ggplot2::geom_tile(ggplot2::aes(fill = value), colour = "white") +
          ggplot2::scale_fill_gradientn(colours = colours) +
          ggplot2::theme(plot.margin = grid::unit(c(0, 0, 0, 0), "lines"))

    plot.list <- list(col = col.plot, row = row.plot, centre = gg)

    invisible(plot.list)
}

#' Show heatmap
#'
#' Assemble and display the components of a cluster heatmap
#'
#' @param plot.list List of plot components returned by makeHeatmap
#' @param col.width Width of the column cluster dendrogram
#' @param row.width Width of the row cluster dendrogram
#'
#' @export
showHeatmap <- function(plot.list, col.width = 0.2, row.width = 0.2) {

    grid::grid.newpage()

    top.layout <- grid::grid.layout(
                    nrow    = 2,
                    ncol    = 2,
                    widths  = grid::unit(c(1 - row.width, row.width), "null"),
                    heights = grid::unit(c(col.width, 1 - col.width), "null"))

    grid::pushViewport(grid::viewport(layout = top.layout))

    if (col.width > 0) {
        print(plot.list$col,
              vp = grid::viewport(layout.pos.col = 1, layout.pos.row = 1))
    }

    if (row.width > 0) {
        print(plot.list$row,
              vp = grid::viewport(layout.pos.col = 2, layout.pos.row = 2))
    }

    print(plot.list$centre +
          ggplot2::theme(axis.line        = ggplot2::element_blank(),
                         axis.text.x      = ggplot2::element_blank(),
                         axis.text.y      = ggplot2::element_blank(),
                         axis.ticks       = ggplot2::element_blank(),
                         axis.title.x     = ggplot2::element_blank(),
                         axis.title.y     = ggplot2::element_blank(),
                         legend.position  = "none",
                         panel.background = ggplot2::element_blank(),
                         panel.border     = ggplot2::element_blank(),
                         panel.grid.major = ggplot2::element_blank(),
                         panel.grid.minor = ggplot2::element_blank(),
                         plot.background  = ggplot2::element_blank()),
          vp = grid::viewport(layout.pos.col = 1, layout.pos.row = 2))

    legend <- heatmapLegend(plot.list$centre)
    grid::pushViewport(grid::viewport(layout.pos.col = 2, layout.pos.row = 1))
    grid::grid.draw(legend)
    grid::upViewport(0)
}

#' Heatmap legend
#'
#' Extract the legend from a ggplot2 heatmap
#'
#' @param heatmap Plot to extract the legend from
#'
#' @return The heatmap legend
heatmapLegend <- function(heatmap) {
    elems  <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(heatmap))
    name   <- which(sapply(elems$grobs, function(x) x$name) == "guide-box")
    legend <- elems$grobs[[name]]
    return(legend)
}
