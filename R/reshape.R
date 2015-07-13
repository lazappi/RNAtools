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

