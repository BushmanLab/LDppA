##' Combine duplicate rows of a table of counts
##'
##' When a table of ocunts has many duplicate rows, computations are
##' sometimes faster if they are combined and a vector of case weights
##' is created.
##' @title wttab Combine Duplicate Rows
##' @param tab matrix or table of counts.  Rows correspond to lineages
##'     and columns to cell types.
##' @return \code{list} with elements \code{tab} - the rows of
##'     \code{unique(tab)}, \code{n} - a vector of row counts, and
##'     \code{data.index} - a mapping of the rows to the \code{tab}
##'     argument to those of \code{wttab(tab)[["tab"]]}.
##' @export
##' @author Charles Berry
wttab <- function(tab){
    tab <- unclass(as.matrix(tab))
    utab <- unique(tab)
    tab.index <- match(
        do.call(paste,as.data.frame(tab)),
        do.call(paste,as.data.frame(utab)))
    list(tab=utab,n=as.vector(table(tab.index)),
         data.index = tab.index)}
