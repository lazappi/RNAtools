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
            dplyr::mutate(row = rownames(data)) %>%
            tidyr::gather(key = row, "value", -row) %>%
            magrittr::set_colnames(c("row", "col", "value"))

    return(long)
}


#' Combine Matrices
#'
#' Combine a list of matrices with the same number of columns into a single long
#' format data.frame
#'
#' @param data.list List of matrices to combine
#'
#' @return Long format data.frame of combined matrices
#'
#' @export
combineMatrices <- function(data.list) {

    long.list <- lapply(data.list, lengthenMatrix)

    for(name in names(long.list)) {
        long.list[[name]]$matrix <- rep(name, nrow(long.list[[name]]))
    }

    combined <- do.call(rbind, long.list)

    return(combined)
}

