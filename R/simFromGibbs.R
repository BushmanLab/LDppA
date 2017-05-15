##' Use the original data and a Gibbs sample to generate a new sample
##' akin to the data.
##'
##' The distribution of cell counts observed for a given composition
##' in a Gibbs sample is used to determine the distribution of cells
##' sampled but not necessarily counted by the cell sorter.  To carry
##' this out, the (zero truncated) distribution of counts of cells for
##' ISs sampled is required.
##'
##' The distribution of sampled cell counts is developed using an EM
##' algorithm.  The distribution of cells observed given the number
##' sampled is binomial with a parameter determined by the
##' composition, the confusion matrix, and the subsampling
##' proportions. 
##'
##' Once the distribution is obtained, a sample of cells is drawn,
##' cells are assigned to cell types, they are \sQuote{sorted}, and
##' subsampled.  The distribution of the number of cells is assumed to
##' be a draw from the symmetric Dirichlet with unit vector as the parameter.
##'
##' The cell counts are rendered in the form given by \code{\link{wttab}}
##'
##' @title Mimic Data From Gibbs Sample
##' @param gibbs The result of \code{\link{estimateComps}}, say.
##' @param params A list such as used by \code{\link{estimateComps}}.
##' @param iterate Which iterate of the Gibbs sampler to use
##' @param wtab The result of \code{\link{wttab}()}
##' @param simpars A list of \code{nr} the largest value to use as
##'     lower abundance, \code{nc} the largest number of observed
##'     cells to include among the lower abundance ISs, and
##'     \code{niter} the number of iterations to use in the EM
##'     fitting.
##' @return A table of counts in the format of \code{wttab()}
##' @importFrom stats dbinom pbinom model.matrix qnbinom
##' @export
##' @author Charles Berry
simFromGibbs <-
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
	    estimateComps(V, eta, alpha, params,
                          wtab, nreps = 1L,
                          nburn = 0L, nthin = 1L,
                          save.last = TRUE,
                          fix.eta=TRUE)[["last"]]
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
    pvec <- eta%*%params[["omega"]]%*%params[["psi"]]

    nonzero.eta <- which(rowSums(zcells)!=0)
    samps <- lapply(nonzero.eta,
		    function(i) samp.Mn(i,zcells,wp.uniq,pvec))

    tab <- simW(samps,lapply(lengths(samps),seq_len),
                eta[nonzero.eta,],params[["omega"]],params[["psi"]])
    wttab(tab)
}

## sample M_n

## note that conditioning on the observed cell counts allows use of the
## negative binomial in place of binomial (conditioning on total counts)

##' 
samp.Mn <- function(ieta,zcells,wp.uniq,pvec,n.iter=10L){
    p.obs <- pvec[ieta]
    zvec <- zcells[ieta,]
    wpvec <- wp.uniq[zvec>0]
    zvec <- zvec[zvec>0]
    ## truncate: do not bother with low probability events
    max.i <- qnbinom(0.9999,max(wpvec)+1,p.obs)+max(wpvec)
    pr.i <- rep(1/max.i,max.i)
    adj.0 <- (1-p.obs)^seq_along(pr.i)/(1-(1-p.obs)^seq_along(pr.i))
    for (iter in seq_len(n.iter)){
        h.i0 <- rep(1,max.i)
        for (i in seq_along(wpvec))
	    h.i0 <-
                h.i0+
                (1+adj.0)* zvec[i]*
                prop.table(pr.i*dnbinom(seq_along(pr.i)-wpvec[i],
                                        wpvec[i]+1,p.obs))
        pr.i <- h.i0/sum(h.i0)
    }
    ## got pr.i
    samps <- rowSums(
        sapply(seq_along(wpvec),
	       function(iwp){
                   rmultinom(1, zvec[iwp],
			     pr.i*
			     dnbinom(seq_along(pr.i)-wpvec[iwp],
				     wpvec[iwp]+1,p.obs))
	       }))
    samp.0 <- rnbinom(max.i,h.i0/(1+adj.0),1-(1-p.obs)^seq_len(max.i))
    samp.0+samps
}
