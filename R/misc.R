finiteMean <- function(x) {
    return(mean(x[is.finite(x)]))
}

finiteVar <- function(x) {
    return(var(x[is.finite(x)]))
}
