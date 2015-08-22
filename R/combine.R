#' Gene Venn Diagram
#'
#' Produce a Venn diagram of significant genes from a list of differential
#' expression results
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

listIntersect <- function(set.list) {

    if (length(set.list) == 0) {
        union.set <- c()
    } else if (length(set.list) == 1) {
        inter.set <- unique(unlist(set.list))
    } else if (length(set.list) == 2) {
        inter.set <- intersect(set.list[[1]], set.list[[2]])
    } else {
        inter.set <- intersect(set.list[[1]], listIntersect(set.list[-1]))
    }

    return(inter.set)
}

listUnion <- function(set.list) {

    if (length(set.list) == 0) {
        union.set <- c()
    } else if (length(set.list) == 1) {
        union.set <- unique(unlist(set.list))
    } else if (length(set.list) == 2) {
        union.set <- union(set.list[[1]], set.list[[2]])
    } else if (length(set.list) > 2) {
        union.set <- union(set.list[[1]], listUnion(set.list[-1]))
    }

    return(union.set)
}

listSetdiff <- function(set.list1, set.list2) {

    set.list1.inter <- listIntersect(set.list1)
    set.list2.union <- listUnion(set.list2)

    diff.set <- setdiff(set.list1.inter, set.list2.union)

    return(diff.set)
}

vennSets <- function(set.list) {

    set.names <- names(set.list)

    combos <- lapply(1:length(set.list),
                     function(j) combn(names(set.list), j, simplify = FALSE))

    combos <- unlist(combos, recursive = FALSE)

    names(combos) <- sapply(combos, function(i) paste0(i, collapse = "-"))

    venn.sets <- lapply(combos,
                        function(i) listSetdiff(set.list[i],
                                               set.list[setdiff(set.names, i)]))

    return(venn.sets)
}

#' Venn Genes
#'
#' Get list of genes for each region in a Venn diagram from a list of
#' differential expression results.
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
                        regulariseResults(name) %>%
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

    regular.data.list <- list()
    gene.lists        <- list()

    for (name in names(data.list)) {

        regular.data <- data.list[[name]] %>%
                        regulariseResults(name)

        gene.list    <- regular.data %>%
                        dplyr::filter(Significance <= alpha) %>%
                        data.frame() %>%
                        magrittr::extract(, "Gene")

        regular.data.list[[name]] <- regular.data
        gene.lists[[name]]        <- gene.list
    }

    summary <- combineMatrices(regular.data.list, lengthen = FALSE) %>%
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


