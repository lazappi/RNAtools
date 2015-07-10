#' Transform counts
#'
#' Apply transformations to raw counts for applications such as plotting
#'
#' @param data    Matrix of raw counts
#' @param methods Vector of transformation methods to apply, choices are "log",
#'                "vst", "rlog" and "cpm".
#'
#' @return List of transformed matrices
#'
#' @export
transformCounts <- function(data, methods) {

    allowed.methods <- c("log", "vst", "rlog", "cpm")

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
                print("Log")
            },

            vst = {
                print("VST")
            },

            rlog = {
                print("rlog")
            },

            cpm = {
                print("CPM")
            },

            stop(paste("Method", method, "not recognised. Allowed methods are:",
                       paste(allowed.methods, collapse = ", ")))
        )
    }
}
