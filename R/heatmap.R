#' Heatmap Dendrogram
#'
#' Produce a dendrogram for a heatmap using ggplot2 and ggdendro
#'
#' @param data Dendrogram data to plot
#' @param rows Boolean, whether the plot is for rows
#'
#' @return ggplot2 object containing the dendrogram
heatmapDendro <- function(data, rows = TRUE) {

    gg <- ggplot2::ggplot() +
          ggplot2::geom_segment(data = ggdendro::segment(data),
                                ggplot2::aes(x = x, y = y,
                                             xend = xend, yend = yend)) +
          ggplot2::labs(x = NULL, y = NULL) +
          ggdendro::theme_dendro()

    if (rows) {
        gg <- gg +
              ggplot2::scale_x_continuous(
                                  expand = c(0.5 / length(data$labels$x), 0)) +
              ggplot2::coord_flip()
    } else {
        gg <- gg +
              ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                                 hjust = 1))
    }

    return(gg)
}

#' Make Heatmap
#'
#' Produce the group of plots required for a heatmap with cluster dendrograms
#' from a matrix
#'
#' @param data Matrix of data to plot
#'
#' @return list of components of the cluster heatmap
#'
#' @export
makeHeatmap <- function(data) {

    colours <- colorRampPalette(RColorBrewer::brewer.pal(9, "GnBu"))(16)

    row.hc <- hclust(dist(data),    "ward.D")
    col.hc <- hclust(dist(t(data)), "ward.D")

    row.dendro <- ggdendro::dendro_data.dendrogram(as.dendrogram(row.hc),
                                                   type = "rectangle")
    col.dendro <- ggdendro::dendro_data.dendrogram(as.dendrogram(col.hc),
                                                   type = "rectangle")

    row.plot <- heatmapDendro(row.dendro, rows = TRUE, labels = FALSE) +
                ggplot2::theme(plot.margin = grid::unit(c(0, 0, 0, 0), "lines"))
    col.plot <- heatmapDendro(col.dendro, rows = FALSE, labels = TRUE) +
                ggplot2::scale_x_continuous(breaks = 1:ncol(data),
                                            labels = colnames(data)) +
                ggplot2::theme(plot.margin = grid::unit(c(0, 0, 0, 0), "lines"))

    row.ord <- match(row.dendro$labels$label, rownames(data))
    col.ord <- match(col.dendro$labels$label, colnames(data))

    data.ord <- data[row.ord, col.ord]
    dimnames(data.ord) <- NULL
    data.ord <- reshape::melt.array(data.ord)

    gg <- ggplot2::ggplot(data.ord, ggplot2::aes(X2, X1)) +
          ggplot2::geom_tile(ggplot2::aes(fill = value), colour = "white") +
          ggplot2::scale_fill_gradientn(colours = colours) +
          ggplot2::labs(x = NULL, y = NULL) +
          ggplot2::scale_x_continuous(expand = c(0, 0)) +
          ggplot2::scale_y_continuous(expand = c(0, 0), breaks = NULL) +
          ggplot2::theme(plot.margin = grid::unit(c(0, 0, 0, 0), "lines"))

    return.list <- list(col = col.plot, row = row.plot, centre = gg)

    invisible(return.list)
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
    legend <- tmp$grobs[[name]]
    return(legend)
}
