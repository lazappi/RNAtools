#' List Normalise
#'
#' Apply the appropriate normalisation to a list of differential expression
#' objects
#'
#' @param data.list List of differential expression objects
#'
#' @return List of normalised objects
#'
#' @export
listNorm <- function(data.list) {

    normalised <- list()

    for (name in names(data.list)) {

        data <- data.list[[name]]

        switch(
            name,

            "edgeR" = {
                norm <- edgeRNorm(data)
            },

            "DESeq" = {
                norm <- deseqNorm(data)
            },

            "DESeq2" = {
                norm <- deseq2Norm(data)
            },

            "voom" = {
                norm <- voomNorm(data)
            },
        )

        normalised[[name]] <- norm
    }

    return(normalised)
}

#' edgeR Normalise
#'
#' Normalise a DGEList object using edgeR
#'
#' @param dge DGEList object to normalise
#'
#' @return DGEList object with normalistion factors
#'
#' @export
edgeRNorm <- function(dge) {

    dge <- edgeR::calcNormFactors(dge)
    dge <- edgeR::estimateCommonDisp(dge)
    dge <- edgeR::estimateTagwiseDisp(dge)

    return(dge)
}


#' DESeq Normalise
#'
#' Normalise a CountDataSet object using DESeq
#'
#' @param count.data CountDataSet object to normalise
#'
#' @return Normalised CountDataSet object
#'
#' @export
deseqNorm <- function(count.data) {

    count.data <- DESeq::estimateSizeFactors(count.data)
    count.data <- DESeq::estimateDispersions(count.data)

    return(count.data)
}


#' DESeq2 Normalise
#'
#' Normalise a DESeqDataSet using DESeq2
#'
#' @param count.data DESeqDataSet object to normalise
#'
#' @return Normalised DESeqDataSet object
#'
#' @export
deseq2Norm <- function(count.data) {

    count.data <- DESeq2::DESeq(count.data)

    return(count.data)
}

#' voom Normalise
#'
#' Normalise a DGEList object using limma-voom
#'
#' @param dge DGEList object to normalise
#'
#' @return Normalised EList object
#'
#' @export
voomNorm <- function(dge) {

    groups <- dge$samples$group
    design <- model.matrix(~ groups)

    dge       <- edgeR::calcNormFactors(dge)
    voom.data <- limma::voom(dge, design)

    return(voom.data)
}
