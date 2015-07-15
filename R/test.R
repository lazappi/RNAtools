#' List Test
#'
#' Apply the appropriate testing procedure to a list of normalised differential
#' expression objects
#'
#' @param data.list List of normalised differential expression objects
#'
#' @return List of test result objects
#'
#' @export
listTest <- function(data.list) {

    tested <- list()

    return(tested)
}

#' edgeR Test
#'
#' Test a normalised DGEList object using edgeR
#'
#' @param dge Normalised DGEList object to test
#'
#' @return
#'
#' @export
edgeRTest <- function(dge) {


    return( )
}


#' DESeq Test
#'
#' Testa a normalised CountDataSet object using DESeq
#'
#' @param count.data Normalised CountDataSet object to test
#'
#' @return
#'
#' @export
deseqTest <- function(count.data) {

    return( )
}


#' DESeq2 Test
#'
#' Test a normalised DESeqDataSet using DESeq2
#'
#' @param count.data Normalised DESeqDataSet object to normalise
#'
#' @return
#'
#' @export
deseq2Test <- function(count.data) {

    return()
}

#' voom Test
#'
#' Test a normalised EList object using limma-voom
#'
#' @param voom.data EList object to test
#'
#' @return
#'
#' @export
voomNorm <- function(voom.data) {


    return( )
}
