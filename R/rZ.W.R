##' Sample Compositions Conditionally on Data
##'
##' The compositions are given by \code{eta} with marginal
##' probabilities of \code{dZ.V(V)}. Samples are drawn conditionally
##' on the cell type counts given by data, \code{wtab}.
##'
##' If \code{gibbs.obj} is given, the result depends entirely on
##' \code{gibbs.obj$last} and all other arguments are ignored.
##' @title rZ.W
##' @param gibbs.obj the value of \code{\link{estimateComps}(...,
##'     save.last=TRUE)}
##' @param wtab (optional) an object such as \code{\link{wttab}}
##'     creates
##' @param eta a matrix of compositions
##' @param V a vector of stick breaking probabilities
##' @param params a list including elements \dQuote{epsilon},
##'     \dQuote{omega}, and \dQuote{psi}.  See
##'     \code{\link{estimateComps}} for details.
##' @export
##' @return a matrix with \code{nrow(eta)} rows and
##'     \code{nrow(wttab$tab)} columns
##' @author Charles Berry
rZ.W <- function(gibbs.obj=NULL, wtab, eta, V, params)
{
    if (is.null(gibbs.obj))
    {
        ## params
        omega <- params$omega
        psi <- params$psi
        eps <- params$epsilon
        K <- ncol(wtab$tab)
        T <- nrow(eta)
        W <- wtab$tab
        Wplus <- rowSums(W)
        Nw <- length(Wplus)
        n <- wtab$n
        ndat <- as.integer(nrow(W))
        omDPsi <- omega%*%diag(psi)

        call.list <-
            list("zmat",
                 prw = double(T),
                 pz = double(T),
                 eps = as.double( eps ),
                 V = as.double(V),
                 eta = as.double( eta ),
                 omdp = as.double(omDPsi),
                 w = as.integer( W ),
                 wp = as.integer( Wplus ),
                 n = as.integer( n ),
                 T = as.integer( T ),
                 ka = as.integer( K ),
                 ko = as.integer( K ),
                 ndat = as.integer( ndat ),
                 zy = integer(T*ndat),
                 etaomdp = double(T*K),
                 workT = double(T))
    } else {
        res.last <- gibbs.obj[["last"]]
        if (is.null(res.last))
            stop("gibbs.obj must have a `last' component")
        call.list <-
            c( list( "zmat" ),
              res.last[c("prw", "pz", "eps", "V", "eta", "omdp", "w",
                         "wp", "n", "T", "ka", "ko", "ndat", "zy",
                         "etaomdp", "workT" )])
    }

    res <- do.call(.C,call.list)

    matrix(res$zy,nrow=res$T)
}
