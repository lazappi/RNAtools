#' Venn Diagram
#'
#' Produce a Venn diagram from a list of differential expression results
#'
#' @param data.list List of results to combine
#' @param alpha     Significance level for selecting genes
#'
#' @return Venn diagram object
#'
#' @importFrom magrittr "%>%"
#'
#' @export
vennDiagram <- function(data.list, alpha = 0.05) {

    gene.lists <- list()

    for (name in names(data.list)) {

        genes <- data.list[[name]] %>%
                 regulariseResults(name) %>%
                 dplyr::filter(Significance <= 0.05)

        gene.lists[[name]] <- genes$Gene
    }

    colours <- RColorBrewer::brewer.pal(n = length(gene.lists), name = "Dark2")

    venn <- VennDiagram::venn.diagram(gene.lists,
                                      filename   = NULL,
                                     height     = 3000,
                                     width      = 3000,
                                     resolution = 500,
                                     col        = "transparent",
                                     fill       = colours,
                                     alpha      = 0.4,
                                     cex        = 1.5,
                                     fontfamily = "sans",
                                     fontface   = "bold",
                                     cat.col    = colours,
                                     cat.cex    = 1.5,
                                     cat.pos    = 0,
                                     cat.dist   = 0.05,
                                     margin     = 0.1)

    return(venn)
}

#' Jaccard Table
#'
#' Produce a table of Jaccard Indices from a list of differential expression
#' results
#'
#' @param data.list List of results to combine
#' @param alpha     Significance level for selecting genes
#'
#' @return dataframe containing Jaccard indices
#'
#' @importFrom magrittr "%>%"
#'
#' @export
jaccardTable <- function(data.list, alpha = 0.05) {

    regular.data.list <- list()

    for (name in names(data.list)) {

        data <- data.list[[name]]

        regular.data <- regulariseResults(data, name)

        regular.data.list[[name]] <- regular.data
    }

    ndata <- length(regular.data.list)

    jaccard.mat <- matrix(nrow = ndata, ncol = ndata)

    n.members <- rep(NA, ndata)

    for(i in 1:ndata) {

        data1 <- regular.data.list[[i]] %>%
                 dplyr::filter(Significance <= alpha)

        data1 <- data1$Gene

        n.members[i] <- length(data1)

        for(j in 1:ndata) {

            data2 <- regular.data.list[[j]] %>%
                     dplyr::filter(Significance <= alpha)

            data2 <- data2$Gene

            jaccard <- length(intersect(data1, data2)) /
                       length(union(data1, data2))
            jaccard.mat[i, j] <- jaccard
        }
    }

    colnames(jaccard.mat) <- names(data.list)
    rownames(jaccard.mat) <- names(data.list)

    jaccard.table <- data.frame(jaccard.mat, Members = n.members)

    return(jaccard.table)
}
