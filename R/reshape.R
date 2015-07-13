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
