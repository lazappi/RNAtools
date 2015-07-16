#' List Test
#'
#' Apply the appropriate testing procedure to a list of normalised differential
#' expression objects
#'
#' @param data.list List of normalised differential expression objects
#' @param group1 First group to test, the reference or control
#' @param group2 Second group to test, the treatment
#'
#' @return List of test result objects
#'
#' @export
listTest <- function(data.list, group1, group2) {

    tested <- list()

    for (name in names(data.list)) {

        data <- data.list[[name]]

        switch(
            name,

            "edgeR" = {
                test <- edgeRTest(data, group1, group2)
            },

            "DESeq" = {
                test <- deseqTest(data, group1, group2)
            },

            "DESeq2" = {
                test <- deseq2Test(data, group1, group2)
            },

            "voom" = {
                test <- voomTest(data, group1, group2)
            },
        )

        tested[[name]] <- test
    }

    return(tested)
}

#' edgeR Test
#'
#' Test a normalised DGEList object using edgeR
#'
#' @param dge Normalised DGEList object to test
#' @param group1 First group to test, the reference or control
#' @param group2 Second group to test, the treatment
#'
#' @return dataframe containing test results
#'
#' @export
edgeRTest <- function(dge, group1, group2) {

    genes.de <- edgeR::exactTest(dge, pair = c(group1, group2))
    top.tags <- edgeR::topTags(genes.de, n = nrow(dge))

    return(top.tags$table)
}

#' DESeq Test
#'
#' Testa a normalised CountDataSet object using DESeq
#'
#' @param count.data Normalised CountDataSet object to test
#' @param group1 First group to test, the reference or control
#' @param group2 Second group to test, the treatment
#'
#' @return dataframe containing test results
#'
#' @export
deseqTest <- function(count.data, group1, group2) {

    results <- DESeq::nbinomTest(count.data, group1, group2)
    results <- results[order(results$padj), ]

    return(results)
}

#' DESeq2 Test
#'
#' Test a normalised DESeqDataSet using DESeq2
#'
#' @param count.data Normalised DESeqDataSet object to normalise
#' @param group1 First group to test, the reference or control
#' @param group2 Second group to test, the treatment
#'
#' @return DESeqResults object containing test results
#'
#' @export
deseq2Test <- function(count.data, group1, group2) {

    results <- DESeq2::results(count.data, contrast = c("group", group1, group2))
    results <- results[order(results$padj), ]

    return(results)
}

#' voom Test
#'
#' Test a normalised EList object using limma-voom
#'
#' @param voom.data EList object to test
#' @param group1 First group to test, the reference or control
#' @param group2 Second group to test, the treatment
#'
#' @return dataframe containing test results
#'
#' @export
voomTest <- function(voom.data, group1, group2) {

    design <- voom.data$design

    fit <- limma::eBayes(limma::lmFit(voom.data, design))
    top.table <- limma::topTable(fit, coef = paste0("groups", group2), n = Inf,
                                 sort.by = "p")

    return(top.table)
}
