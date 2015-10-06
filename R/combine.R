#' Gene Venn Diagram
#'
#' Produce a Venn diagram of significant genes from a list of regularised
#' differential expression results
#'
#' @param data.list List of results to combine
#' @param alpha     Significance level for selecting genes
#'
#' @return Venn diagram object
#'
#' @importFrom magrittr "%>%"
#'
#' @export
geneVenn <- function(data.list, alpha = 0.05) {

    gene.lists <- list()

    for (name in names(data.list)) {

        genes <- data.list[[name]] %>%
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
#' Produce a table of Jaccard Indices from a list of regularised differential
#' expression results
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

    ndata <- length(data.list)

    jaccard.mat <- matrix(nrow = ndata, ncol = ndata)

    n.members <- rep(NA, ndata)

    for(i in 1:ndata) {

        data1 <- data.list[[i]] %>%
                 dplyr::filter(Significance <= alpha)

        data1 <- data1$Gene

        n.members[i] <- length(data1)

        for(j in 1:ndata) {

            data2 <- data.list[[j]] %>%
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

#' Venn Genes
#'
#' Get list of genes for each region in a Venn diagram from a list of
#' regularised differential expression results.
#'
#' @param data.list List of results to combine
#' @param alpha     Significance level for selecting genes
#'
#' @return List of vectors of genes
#'
#' @importFrom magrittr "%>%"
#'
#' @export
vennGenes <- function(data.list, alpha = 0.05) {

    gene.lists <- list()

    for (name in names(data.list)) {

        regular.data <- data.list[[name]] %>%
                        dplyr::filter(Significance <= alpha)

        gene.lists[[name]] <- regular.data$Gene
    }

    venn.genes <- vennSets(gene.lists)

    return(venn.genes)
}

#' Gene Summary
#'
#' Summarise the results from multiple differential expression methods by gene.
#' Returns a table with summary statistics and a count of the number of methods
#' that found the gene to be differentially expressed at a give significance
#' level.
#'
#' @param data.list List of results to summarise
#' @param alpha     Significance level for selecting genes
#'
#' @return data.table containing summary statistics
#'
#' @importFrom magrittr "%>%"
#'
#' @export
geneSummary <- function(data.list, alpha = 0.05) {

    gene.lists <- list()

    for (name in names(data.list)) {

        data      <- data.list[[name]]

        gene.list <- data %>%
                     dplyr::filter(Significance <= alpha) %>%
                     data.frame() %>%
                     magrittr::extract(, "Gene")

        gene.lists[[name]]        <- gene.list
    }

    summary <- combineMatrices(data.list, lengthen = FALSE) %>%
               dplyr::group_by(Gene) %>%
               dplyr::summarise(meanFC   = mean(FoldChange),
                                varFC    = var(FoldChange),
                                fMeanFC  = finiteMean(FoldChange),
                                fVarFC   = finiteVar(FoldChange),
                                meanPVal = mean(pValue),
                                varPVal  = var(pValue),
                                meanSig  = mean(Significance),
                                varSig   = var(Significance))

    de.count <- gene.lists %>% unlist %>% table
    de.count <- de.count[summary$Gene]
    de.count[is.na(de.count)] <- 0

    summary$DECount <- de.count

    return(summary)

}


