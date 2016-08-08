##' Faster versions of \code{prop.table} and \code{margin.table}
##'
##' Avoid lengthy loops when possible
##' @title Margin Sums/Proportions
##' @param x table
##' @param margin index, or vector of indices to generate margin for
##' @return same as \code{prop.table} or \code{margin.table}
##' @author Charles Berry
quick.margin.table <- function(x, margin = NULL)
{
    if (is.null(margin)) return(sum(x))
    dmx <- dim(x)
    if (prod(dmx[margin])<prod(dmx[-margin])) return(margin.table(x,margin))

    dmx[margin] <- NA
    mgrid <-
        expand.grid(lapply(dmx,
                           function(x) if (is.na(x)) TRUE else 1:x))
    extractx <- function(...) x[...]
    Reduce("+",
           lapply(1:nrow(mgrid),
                  function(y) do.call( extractx, mgrid[y,,drop=FALSE])))
}

quick.prop.table <- function(x, margin=NULL){
    if (length(margin)) 
        sweep(x, margin, quick.margin.table(x, margin), "/", check.margin = FALSE)
    else x/sum(x)
}
