##' Use the original data and a Gibbs sample to generate a new sample
##' akin to the data.
##'
##' The distribution of cell counts observed for a given composition
##' in a Gibbs sample is used to determine the distribution of cells
##' sampled but not necessarily counted by the cell sorter.  To carry
##' this out, the (zero truncated) distribution of counts of cells for
##' ISs sampled is required.
##'
##' For lower abundance ISs, the distribution of sampled cell counts
##' is developed using an EM algorithm.  The distribution of cells
##' observed given the number sampled is binomial with a parameter
##' determined by the composition, the confusion matrix, and the
##' subsampling proportions.  The proportions are smoothed to temper
##' issues with small cell counts.
##'
##' Once the distribution is obtained, a sample of cells is drawn,
##' cells are assigned to cell types, they are \sQuote{sorted}, and
##' subsampled.  For higher abundance ISs, the distribution of the
##' number of cells is assumed to be locally flat, which leads to a
##' negative binomial distribution with the complement of the binomial
##' probability and one more success than the number of observed cells.
##'
##' The cell counts are rendered in the form given by \code{\link{wttab}}
##'
##' @title Mimic Data From Gibbs Sample
##' @param gibbs The result of \code{\link{ldppa.gibbs.3}}, say.
##' @param params A list such as used by \code{\link{ldppa.gibbs.3}}.
##' @param iterate Which iterate of the Gibbs sampler to use
##' @param wtab The result of \code{\link{wttab}()}
##' @param simpars A list of \code{nr} the largest value to use as
##'     lower abundance, \code{nc} the largest number of observed
##'     cells to include among the lower abundance ISs, and
##'     \code{niter} the number of iterations to use in the EM
##'     fitting.
##' @return A table of counts in the format of \code{wttab()}
##' @importFrom stats dbinom pbinom model.matrix
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
    wp.factor <- factor(pmin(last[["wp"]], nc+1L),
                        1:(nc+1L))
    zcells <- zmat%*%unname(model.matrix(~0+wp.factor))
    pvec <- eta%*%params[["omega"]]%*%params[["psi"]]

    ## fit distrn of low abundance sites
    fit.z <-
        function(p.obs,zc,nr,niter=simpars[["niter"]])
    {
        nc <- length(zc)
        mat <-
            t(sapply(1:nr,
                     function(x) dbinom(1:nc,x,p.obs)))
        m0 <- rowSums(mat)
        po <- rep(1,nr)/nr
        emat <- prop.table(diag(po)%*%mat,2)%*%diag(zc)
        etot <- rowSums(emat)/m0
        for (i in 1:niter){
            po <- prop.table(etot+1)
            emat <- prop.table(diag(po)%*%mat,2)%*%diag(zc)
            etot <- rowSums(emat)/m0
        }
        emat
    }

    fits <- lapply(1:nrEta,function(r) fit.z(pvec[r],zcells[r,-1-nc],nr))

    nsampFun <- function(r){
        n <- rowSums(apply(fits[[r]],2,
                           function(x) {
			       colsum <- round(sum(x))
			       if (colsum!=0)
                                   rmultinom(1,colsum,x)
			       else
                                   rep(0,length(x))}))
        m <-
            mapply(function(x,y)
                if (x==0)
                    x
                else
                    rnbinom(1,x,pbinom(0,y,pvec[r],lower.tail=FALSE)),
                n,1:nr)
        m+n
    }

    nsamp.reps <-
        lapply(1:nrEta,nsampFun)                                                


    smaller.tab <-
        simW(nsamp.reps,
             rep(list(1:nr),nrEta),
             eta, params[["omega"]], params[["psi"]])

    ## high abundance sites 
    wp.bigger <- wp.factor == nc+1
    wp.big.value <- last[["wp"]][wp.bigger]
    zmat.big <- zmat[,wp.bigger]
    zmat.counts <- zmat.big[zmat.big!=0]
    zmat.etaRow <- row(zmat.big)[zmat.big!=0]
    zmat.wp <- wp.big.value[col(zmat.big)[zmat.big!=0]]
    bigger <- rnbinom(sum(zmat.counts),
		      1+rep(zmat.wp,zmat.counts),
		      rep(pvec[zmat.etaRow],zmat.counts)) +
        rep(zmat.wp,zmat.counts)
    Nsamps <- split(bigger,factor(rep(zmat.etaRow,zmat.counts),1:nrEta))

    bigger.tab <- simW(lapply(Nsamps,function(x) rep(1,length(x))) ,
		       Nsamps, eta, params[["omega"]], params[["psi"]])
    wttab(rbind(smaller.tab,bigger.tab))
}
