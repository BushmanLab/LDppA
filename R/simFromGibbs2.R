##' Use the original data and a Gibbs sample to generate a new sample
##' akin to the data. Fix the distribution of observed cell counts.
##'
##' The distribution of cell counts observed for a given composition
##' in a Gibbs sample is used to determine the distribution of cells
##' sampled by type.
##'
##' The cell counts are rendered in the form given by \code{\link{wttab}}
##'
##' @title Mimic Data From Gibbs Sample fixing \code{rowSums(wtab[["tab"]])}
##' @param gibbs The result of \code{\link{ldppa.gibbs.3}}, say.
##' @param params A list such as used by \code{\link{ldppa.gibbs.3}}.
##' @param iterate Which iterate of the Gibbs sampler to use
##' @param wtab The result of \code{\link{wttab}()}
##' @param simpars not used
##' @return A table of counts in the format of \code{wttab()}
##' @importFrom stats dbinom pbinom model.matrix qnbinom
##' @export
##' @author Charles Berry
simFromGibbs2 <-
    function(gibbs,params,iterate=NULL,wtab=NULL,
             simpars=list(nr=60L,nc=15L,niter=100L))
{
    niters <- length(gibbs[["monitors"]])
    nr <- simpars[["nr"]]
    nc <- simpars[["nc"]]

    ## default to last iterate
    if (is.null(iterate)) 
        iterate <- niters

    V <- gibbs[["Vs"]][[iterate]]
    eta <- gibbs[["etas"]][[iterate]]
    alpha <- gibbs[["alphas"]][[iterate]]

    ## get result of a single iterate
    last <-
        if (iterate != niters || is.null(gibbs[["last"]]))
            ldppa.gibbs.4(V, eta, alpha, params,
                          wtab, nreps = 1L,
                          nburn = 0L, nthin = 1L,
                          save.last = TRUE)[["last"]]
        else
            gibbs[["last"]]

    nrEta <- last[["T"]]
    nrTab <- last[["ndat"]]
    if (is.null(wtab))
        wtab <- list(tab=matrix(last[["w"]],nrow=nrTab),
                     n=last[["n"]])    
    zmat <- matrix(last[["zy"]],nrow=nrEta)
    wp.uniq <- unique(last[["wp"]])
    wp.factor <- factor(last[["wp"]], wp.uniq)
    zcells <- zmat%*%unname(model.matrix(~0+wp.factor))

    zelts <- which(zcells!=0,arr.ind=TRUE)

    wtt <- lapply(1:nrow(zelts),
                  function(x) {
		      xr <- zelts[x,,drop=FALSE]
		      rmultinom(zcells[xr],
                                wp.uniq[xr[2]],
                                eta[xr[1],]%*%
                                params[["omega"]]%*%
                                diag(params[["psi"]]))})
    ttab <- unlist(wtt,use.names=FALSE)
    dim(ttab) <- c(ncol(eta),length(ttab)/ncol(eta))
    wttab(t(ttab))
}
