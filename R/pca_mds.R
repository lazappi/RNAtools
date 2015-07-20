#' List PCA
#'
#' Produce PCA plots from a list of matrices
#'
#' @param data.list List of matrices to plot
#' @param top       Number of rows with highes variance to select for plotting
#' @param groups    Vector of groups assigned to sample columns
#'
#' @return List of ggplot2 objects containing PCA plots
#'
#' @importFrom magrittr "%>%"
#'
#' @export
listPCA <- function(data.list, top = nrow(data.list[[1]]),
                    groups = colnames(data.list[[1]])) {

    plots <- list()

    # Produce individual plots
    for (name in names(data.list)) {

        gg <- samplePCA(data.list[[name]], top = top , group = groups) +
              ggplot2::ggtitle(paste(name," - PC1 vs PC2, top", top,
                                     "variable genes"))

        plots[[name]] <- gg
    }

    # Produce combined plots
    gg <- data.list %>%
          lapply(samplePCA, top = top, group = groups, plot = FALSE) %>%
          combineMatrices(lengthen = FALSE) %>%
          ggplot2::ggplot(ggplot2::aes(x = PC1, y = PC2, colour = Group,
                                     label = Sample)) +
          ggplot2::geom_text() +
          ggplot2::facet_wrap(~ matrix) +
          ggplot2::ggtitle(paste("PC1 vs PC2, top", top, "variable genes")) +
          ggplot2::xlab("PC1") +
          ggplot2::ylab("PC2")

    plots[["combined"]] <- gg

    return(plots)
}

#' Sample PCA Plot
#'
#' Produce a PCA plot from a matrix using ggplot2
#'
#' @param data   Matrix of data to plot
#' @param top    Number of rows with highest variance to select for plotting
#' @param groups Vector of groups assigned to sample columns
#' @param plot   Boolean, if true return plot, if false return plot data
#'
#' @return ggplot2 object containing the PCA plot, or the dataframe of plot data
#'
#' @importFrom magrittr "%>%"
#'
#' @export
samplePCA <- function(data, top = nrow(data), groups = colnames(data),
                      plot = TRUE) {

    # Select high variance genes
    top.data <- data %>%
                data.frame %>%
                dplyr::mutate(var = genefilter::rowVars(data)) %>%
                dplyr::top_n(top, var) %>%
                dplyr::select(-var)

    # Compute PCA
    PCA.data <- prcomp(t(top.data), scale = FALSE)

    # Calculate variance
    percent.var <- round(100 * PCA.data$sdev ^ 2 / sum(PCA.data$sdev ^ 2), 1)

    # Reshape data from plotting
    plot.data <- data.frame(PC1    = PCA.data$x[, 1],
                            PC2    = PCA.data$x[, 2],
                            Sample = colnames(data),
                            Group  = groups)

    gg <- ggplot2::ggplot(plot.data,
                          ggplot2::aes(x = PC1, y = PC2, colour = Group,
                                       label = Sample)) +
          ggplot2::geom_text() +
          ggplot2::ggtitle(paste("PC1 vs PC2, top", top, "variable genes")) +
          ggplot2::xlab(paste0("PC1, VarExp: ",
                               round(percent.var[1], 4), "%")) +
          ggplot2::ylab(paste0("PC2, VarExp: ",
                               round(percent.var[2], 4), "%"))
    if (plot) {
        return(gg)
    } else {
        return(plot.data)
    }
}

#' Sample MDS plot
#'
#' Produce a MDS plot from a matrix using ggplot2. Based on the limma function.
#'
#' @param data   Matrix of data to plot
#' @param top    Number of rows with highest deviance to select for plotting
#' @param groups Vector of groups assigned to sample columns
#' @param plot   Boolean, if true return plot, if false return plot data
#' @param selection Select rows in a "pariwise" manner between samples or
#'                  "common" across all samples
#'
#' @return ggplot2 object containing the PCA plot, or the dataframe of plot data
#'
#' @importFrom magrittr "%>%"
#'
#' @export
sampleMDS <- function(data, top = nrow(data), groups = colnames(data),
                    plot = TRUE, selection = c("pairwise", "common")) {

    if (missing(selection)) {
        selection <- "pairwise"
    } else {
        selection <- match.arg(selection)
    }

    # Create matrix for storing distances
    dists <- matrix(0, nrow = ncol(data), ncol = ncol(data),
                    dimnames = list(colnames(data), colnames(data)))

    if (selection == "pairwise") {

        top.idx <- ncol(data) - top + 1

        # For each pair of samples
        for (i in 2:ncol(data)) {
            for (j in 1:(i - 1)) {

                # Calculate distance high variance genes for this pair
                sqr.dist <- sort((data[, i] - data[, j]) ^ 2, decreasing = TRUE)
                sqr.dist <- sqr.dist[1:top]

                dists[i, j] <- sqrt(mean(sqr.dist))
            }
        }

        axis.label <- "Leading logFC dim"

    } else {

        # Find high variance genes
        top.data <- data %>%
                    data.frame %>%
                    dplyr::mutate(sqr.dev =
                                  rowMeans((data - rowMeans(data)) ^ 2)) %>%
                    dplyr::top_n(top, sqr.dev) %>%
                    dplyr::select(-sqr.dev)

        # Calculate ditances
        for (i in 2:ncol(data)) {
            dists[i, 1:(i - 1)] <- sqrt(colMeans(
                                            (top.data[, i] -
                                                 top.data[, 1:(i - 1),
                                                          drop = FALSE]) ^ 2))
        }

        axis.label <- "Principal Component"
    }

    # Calculate MDS
    MDS.data <- cmdscale(as.dist(dists), k = 2)

    # Reshape data for plotting
    plot.data <- data.frame(X      = MDS.data[, 1],
                            Y      = MDS.data[, 2],
                            Sample = colnames(data),
                            Group  = groups)

    gg <- ggplot2::ggplot(plot.data,
                          ggplot2::aes(x = X, y = Y, colour = Group,
                                       label = Sample)) +
          ggplot2::geom_text() +
          ggplot2::ggtitle(paste("MDS Plot, top", top, selection,
                                 "variable genes")) +
          ggplot2::xlab(paste(axis.label, 1)) +
          ggplot2::ylab(paste(axis.label, 2))

    if (plot) {
        return(gg)
    } else {
        return(plot.data)
    }
}

#' List MDS
#'
#' Produce MDS plots from a list of matrices
#'
#' @param data.list List of matrices to plot
#' @param top       Number of rows with highes variance to select for plotting
#' @param groups    Vector of groups assigned to sample columns
#' @param selection Select rows in a "pariwise" manner between samples or
#'                  "common" across all samples
#'
#' @return List of ggplot2 objects containing MDS plots
#'
#' @importFrom magrittr "%>%"
#'
#' @export
listMDS <- function(data.list, top = nrow(data.list[[1]]),
                    groups = colnames(data.list[[1]]), selection = "pairwise") {

    plots <- list()

    # Produce indiviual plots
    for (name in names(data.list)) {

        gg <- sampleMDS(data.list[[name]], top = top , group = groups,
                      selection = selection) +
              ggplot2::ggtitle(paste(name," - MDS Plot, top", top,
                                     selection, "variable genes"))

        plots[[name]] <- gg
    }

    # Set appropriate lables
    if (selection == "pairwise") {
        axis.label <- "Leading logFC dim"
    } else {
        axis.label <- "Principal Component"
    }

    # Produce combined plot
    gg <- data.list %>%
        lapply(sampleMDS, top = top, group = groups, plot = FALSE,
               selection = selection) %>%
        combineMatrices(lengthen = FALSE) %>%
        ggplot2::ggplot(ggplot2::aes(x = X, y = Y, colour = Group,
                                     label = Sample)) +
        ggplot2::geom_text() +
        ggplot2::facet_wrap(~ matrix) +
        ggplot2::ggtitle(paste("MDS Plots, top", top, selection,
                               "variable genes")) +
        ggplot2::xlab(paste(axis.label, 1)) +
        ggplot2::ylab(paste(axis.label, 2)) #+
        #ggplot2::theme(axis.title   = ggplot2::element_text(size = 20),
        #               axis.text    = ggplot2::element_text(size = 15),
        #               plot.title   = ggplot2::element_text(size = 30,
        #                                                    face = "bold"),
        #               legend.text  = ggplot2::element_text(size = 15),
        #               legend.title = ggplot2::element_text(size = 15))

    plots[["combined"]] <- gg

    return(plots)
}
