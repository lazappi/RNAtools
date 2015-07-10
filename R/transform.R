#' Transform counts
#'
#' Apply transformations to raw counts for applications such as plotting
#'
#' @param data    Matrix of raw counts
#' @param methods Vector of transformation methods to apply, choices are "log",
#'                "vst", "rlog" and "cpm".
#' @param verbose Boolean, whether to print messages showing progress
#'
#' @return List of transformed matrices
#'
#' @export
transformCounts <- function(data, methods, verbose = TRUE) {

    allowed.methods <- c("log", "vst", "rlog", "cpm", "logCPM")

    if (missing(data)) {
        stop("No data provided")
    } else if (missing(methods)) {
        stop("Methods not provided")
    }

    methods <- unique(methods)

    transformed <- list()

    for (method in methods) {
        switch(
            method,

            log = {
                if (verbose) {message("Applying log transformation...")}
                transformed$log <- logTransform(data)
            },

            vst = {
                if (verbose) {
                    message("Applying variance stabilising transformation...")
                }
                transformed$vst <- vstTransform(data)
            },

            rlog = {
                if (verbose) {
                    message("Applying regularised log transformation...")
                }
                transformed$rlog <- rlogTransform(data)
            },

            cpm = {
                if (verbose) {message("Applying CPM transformation...")}
                transformed$cpm <- cpmTransform(data)
            },

            logCPM = {
                if (verbose) {message("Applying log CPM transformation...")}
                transformed$logCPM <- logCPMTransform(data)
            },

            stop(paste("Method", method, "not recognised. Allowed methods are:",
                       paste(allowed.methods, collapse = ", ")))
        )
    }

    # Sort list by method name
    transformed <- transformed[order(names(transformed))]

    return(transformed)
}

#' Log transform counts
#'
#' Apply the log transformation to raw counts
#'
#' @param data Matrix of raw counts
#'
#' @return Matrix of transformed data
#'
#' @export
logTransform <- function(data) {
    transformed <- log2(data + 1)

    return(transformed)
}

#' VST transform counts
#'
#' Apply DESeq's variance stabilising transformation to raw counts
#'
#' @param data Matrix of raw counts
#'
#' @return Matrix of transformed data
#'
#' @export
vstTransform <- function(data) {

    # Create CountDataSet object with dummy conditions
    count.data <- DESeq::newCountDataSet(data, conditions = rep(1, ncol(data)))
    count.data <- DESeq::estimateSizeFactors(count.data)
    count.data <- DESeq::estimateDispersions(count.data, method = "blind")

    transformed <-DESeq::getVarianceStabilizedData(count.data)

    return(transformed)
}

#' Regularised log transform counts
#'
#' Apply DESeq2's regularised log transformation to raw counts
#'
#' @param data Matrix of raw counts
#'
#' @return Matrix of transformed data
#'
#' @export
rlogTransform <- function(data) {

    transformed <- DESeq2::rlogTransformation(counts)
    rownames(transformed) <- rownames(data)

    return(transformed)
}

#' Counts Per Million transform counts
#'
#' Apply edgeR's Counts Per Million transformation to raw counts
#'
#' @param data Matrix of raw counts
#'
#' @return Matrix of transformed data
#'
#' @export
cpmTransform <- function(data) {

    data <- edgeR::DGEList(data)
    transformed <- edgeR::cpm(data)

    return(transformed)
}

#' Log Counts Per Million transform counts
#'
#' Apply edgeR's Counts Per Million transformation to raw counts, then take the
#' log
#'
#' @param data Matrix of raw counts
#'
#' @return Matrix of transformed data
#'
#' @export
logCPMTransform <- function(data) {

    data <- edgeR::DGEList(data)
    transformed <- edgeR::cpm(data, log = TRUE)

    return(transformed)
}
