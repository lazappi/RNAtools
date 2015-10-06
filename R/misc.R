finiteMean <- function(x) {
    return(mean(x[is.finite(x)]))
}

finiteVar <- function(x) {
    return(var(x[is.finite(x)]))
}

listIntersect <- function(set.list) {

    if (length(set.list) == 0) {
        union.set <- c()
    } else if (length(set.list) == 1) {
        inter.set <- unique(unlist(set.list))
    } else if (length(set.list) == 2) {
        inter.set <- intersect(set.list[[1]], set.list[[2]])
    } else {
        inter.set <- intersect(set.list[[1]], listIntersect(set.list[-1]))
    }

    return(inter.set)
}

listUnion <- function(set.list) {

    if (length(set.list) == 0) {
        union.set <- c()
    } else if (length(set.list) == 1) {
        union.set <- unique(unlist(set.list))
    } else if (length(set.list) == 2) {
        union.set <- union(set.list[[1]], set.list[[2]])
    } else if (length(set.list) > 2) {
        union.set <- union(set.list[[1]], listUnion(set.list[-1]))
    }

    return(union.set)
}

listSetdiff <- function(set.list1, set.list2) {

    set.list1.inter <- listIntersect(set.list1)
    set.list2.union <- listUnion(set.list2)

    diff.set <- setdiff(set.list1.inter, set.list2.union)

    return(diff.set)
}

vennSets <- function(set.list) {

    set.names <- names(set.list)

    combos <- lapply(1:length(set.list),
                     function(j) combn(names(set.list), j, simplify = FALSE))

    combos <- unlist(combos, recursive = FALSE)

    names(combos) <- sapply(combos, function(i) paste0(i, collapse = "-"))

    venn.sets <- lapply(combos,
                        function(i) listSetdiff(set.list[i],
                                                set.list[setdiff(set.names, i)]))

    return(venn.sets)
}
