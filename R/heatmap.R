

dendroPlot <- function(data, rows = TRUE, labels = TRUE) {

    # y.range <- range(data$segments$y)
    # y.diff  <- y.range[2] - y.range[1]
    # n.char  <- max(nchar(as.character(data$labels$label)))
    # angle   <- if (rows) {0} else {90}

    gg <- ggplot() +
          geom_segment(data = segment(data),
                       aes(x = x, y = y, xend = xend, yend = yend)) +
          labs(x = NULL, y = NULL) +
          theme_dendro()

    if (rows) {
        gg <- gg +
              scale_x_continuous(expand = c(0.5 / length(data$labels$x), 0)) +
              coord_flip()
    } else {
        gg <- gg + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    }

    return(gg)
}

ggheatmap <- function(data) {

    colours <- colorRampPalette(brewer.pal(9, "GnBu"))(16)

    row.hc <- hclust(dist(x), "ward")
    col.hc <- hclust(dist(x), "ward")

    row.dendro <- dendro_data(as.dendrogram(row.hc), type = "rectangle")
    col.dendro <- dendro_data(as.dendrogram(col.hc), type = "rectangle")

    row.plot <- dendroPlot(row.dendro, rows = TRUE, labels = FALSE) +
                theme(plot.margin = unit(c(0, 0, 0, 0), "lines"))
    col.plot <- dendroPlot(col.dendro, rows = FALSE, labels = TRUE) +
                scale_x_continuous(breaks = 1:ncol(data),
                                   labels = colnames(data)) +
                theme(plot.margin = unit(c(0, 0, 0, 0), "lines"))

    row.ord <- match(row.dendro$labels$label, rownames(data))
    col.ord <- match(col.dendro$labels$label, colnames(data))

    data.ord <- data[row.ord, col.ord]
    dimnames(data.ord) <- NULL
    data.ord <- melt(data.ord)

    gg <- ggplot(data.ord, aes(X2, X1)) +
          geom_tile(aes(fill = value), colour = "white") +
          scale_till_gradientn(colours = colours) +
          labs(x = NULL, y = NULL) +
          scale_x_continuous(expand = c(0, 0)) +
          scale_y_continuous(expand = c(0, 0), breaks = NULL) +
          theme(plot.margin = unit(c(0, 0, 0, 0), "lines"))

    return.list <- list(col = col.plot, row = row.plot, centre = gg)

    invisible(return.list)
}

showHeatmap <- function(plot.list, col.width = 0.2, row.width = 0.2) {

    gt <- gtable(widths  = unit(c(1 - row.width, row.width), "null"),
                 heights = unit(c(1 - col.width, col.width), "null"))

    grobs <- lapply(plot.list, function(p) gg(gg(p)))

    if (col.width > 0) {
        gt <- gtable_add_grob(gt, grobs$col[[]], 1, 1)
    }

    print(plot.list$col, vp = viewport(layout.pos.col = 1, layout.pos.row = 1))

    if (row.width > 0) {
        print(plot.list$row,
              vp = viewport(layout.pos.col = 2, layout.pos.row = 2))
    }

    print(plot.data$centre +
          theme(axis.line = element_blank(),
                axis.text.x = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks = element_blank(),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                legend.position = "none",
                panel.background = element_blank(),
                panel.border = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                plot.background = element_blank()),
          vp = viewport(layout.pos.col = 1, layout.pos.row = 2))

    legend <- g_legend(plot.list$centre)
    pushViewport(viewport(layout.pos.col = 2, layout.pos.row = 1))
    grid.draw(legend)
    upViewport(0)
}
