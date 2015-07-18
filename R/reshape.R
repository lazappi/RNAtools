#' Lengthen
#'
#' Convert a matrix with informative rowanames to a long format data.frame
#'
#' @param data Matrix to convert to long format
#'
#' @return Long format data.frame
#'
#' @importFrom magrittr "%>%"
#'
#' @export
lengthenMatrix <- function(data) {

    long <- data %>%
            data.frame %>%
            dplyr::add_rownames() %>%
            tidyr::gather(key = rowname, "value", -rowname) %>%
            magrittr::set_colnames(c("row", "col", "value"))

    return(long)
}

#' Combine Matrices
#'
#' Combine a list of matrices with the same number of columns into a single long
#' format data.frame
#'
#' @param data.list List of matrices to combine
#' @param lengthen  Boolean, whether to apply lengthenMatrix before combining
#'
#' @return Long format data.frame of combined matrices
#'
#' @export
combineMatrices <- function(data.list, lengthen = TRUE) {

    if (lengthen) {
        data.list <- lapply(data.list, lengthenMatrix)
    }

    for(name in names(data.list)) {
        data.list[[name]]$matrix <- rep(name, nrow(data.list[[name]]))
    }

    combined <- do.call(rbind, data.list)

    return(combined)
}

#' Counts to Objects
#'
#' Convert a count matrix to objects required by various differential expression
#' testing packages
#'
#' @param data    Matrix of counts
#' @param groups  Vector of groups for each sample column
#' @param filter  Boolean, whether to use HTSFilter for voom
#' @param methods Vector of object types to convert to
#' @param verbose Boolean, whether to print messages showing progress
#'
#' @return List of converted objects
#'
#' @export
counts2Objects <- function(data, groups, filter,
                           methods = c("edgeR", "DESeq", "DESeq2", "voom"),
                           verbose = TRUE) {

    methods <- unique(methods)

    objects <- list()

    for (method in methods) {
        switch(
            method,

            edgeR = {
                if (verbose) {message("Creating object for edgeR...")}
                objects$edgeR <- counts2edgeR(data, groups)
            },

            DESeq = {
                if (verbose) message("Creating object for DESeq...")
                objects$DESeq <- counts2DESeq(data, groups)
            },

            DESeq2 = {
                if (verbose) message("Creating DESeq2 object...")
                objects$DESeq2 <- counts2DESeq2(data, groups)
            },

            voom = {
                if (verbose) {message("Creating voom object...")}
                objects$voom <- counts2voom(data, groups, filter)
            },

            stop(paste("Method", method, "not recognised. Allowed methods are:",
                       "edgeR/DESeq/DESeq2/voom"))
        )
    }

    # Sort list by method name
    objects <- objects[order(names(objects))]

    return(objects)
}

#' Counts to edgeR
#'
#' Convert a count matrix to an edgeR DGEList object
#'
#' @param data    Matrix of counts
#' @param groups  Vector of groups for each sample column
#'
#' @return DGEList object containing counts
#'
#' @export
counts2edgeR <- function(data, groups) {

    dge <- edgeR::DGEList(data, group = groups)

    return(dge)
}

#' Counts to voom
#'
#' Convert a count matrix to an edgeR DGEList object for voom
#'
#' @param data   Matrix of counts
#' @param groups Vector of groups for each sample column
#' @param filter Boolean, whether to apply HTSFilter
#'
#' @return DGEList object containing counts
#'
#' @export
counts2voom <- function(data, groups, filter) {

    dge <- edgeR::DGEList(data, group = groups)

    if (filter) {
        filtered <- HTSFilter::HTSFilter(dge, plot = FALSE)
        dge <- filtered$filteredData

        message(paste("HTSFilter threshold:", filtered$s,
                      "Genes Filtered:", nrow(data) - nrow(dge$counts)))

    }

    return(dge)
}

#' Counts to DESeq
#'
#' Convert a count matrix to a DESeq CountDataSeq object
#'
#' @param data    Matrix of counts
#' @param groups  Vector of groups for each sample column
#'
#' @return CountDataSet object containing counts
#'
#' @export
counts2DESeq <- function(data, groups) {

    count.data <- DESeq::newCountDataSet(data, conditions = groups)

    return(count.data)
}

#' Counts to DESeq2
#'
#' Convert a count matrix to a DESeq2 DESeqDataSet object
#'
#' @param data    Matrix of counts
#' @param groups  Vector of groups for each sample column
#'
#' @return CountDataSet object containing counts
#'
#' @export
counts2DESeq2 <- function(data, groups) {

    groups <- data.frame(group = groups)

    count.data <- DESeq2::DESeqDataSetFromMatrix(data,
                                                 colData = groups,
                                                 design = ~ group)

    return(count.data)
}


regulariseResults <- function(results,
                              method = c("edgeR", "DESeq", "DESeq2", "voom") ) {

    # Check that a valid method has been given
    if (missing(method)) {
        stop("Differential expression method must be specified")
    } else {
        method <- match.arg(method)
    }

    switch(
        method,

        edgeR = {
            regular <- results %>%
                       dplyr::add_rownames() %>%
                       dplyr::select(Gene         = rowname,
                                     FoldChange   = logFC,
                                     Abundance    = logCPM,
                                     Significance = FDR)
        },

        DESeq = {
            regular <- results %>%
                       dplyr::select(Gene         = id,
                                     FoldChange   = log2FoldChange,
                                     Abundance    = baseMean,
                                     Significance = padj) %>%
                       dplyr::mutate(Abundance = log2(Abundance))
        },

        DESeq2 = {
            regular <- results %>%
                       data.frame %>%
                       dplyr::add_rownames() %>%
                       dplyr::select(Gene         = rowname,
                                     FoldChange   = log2FoldChange,
                                     Abundance    = baseMean,
                                     Significance = padj) %>%
                       dplyr::mutate(Abundance = log2(Abundance))
        },

        voom = {
            regular <- results %>%
                       dplyr::add_rownames() %>%
                       dplyr::select(Gene         = rowname,
                                     FoldChange   = logFC,
                                     Abundance    = AveExpr,
                                     Significance = adj.P.Val)
        }
    )

    return(regular)
}
