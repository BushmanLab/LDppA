##' Estimate Cell Type Composition Weights Given Eta
##'
##' This maximizes a penalized likelihood using the EM algorithm
##' @title fixed Eta EM
##' @param wtab object as created by \code{\link{wttab}}
##' @param eta a matrix whose rows are compositional probabilities
##' @param params list with matrix omega and vector psi
##' @param alpha positive numeric atom. Use \code{alpha==0.0} for
##'     maximum likelihood (which forces \code{weights==1}
##' @param p.init numeric vector of mixing proportions
##' @param weights non0negative number vector with one element for
##'     each composition.  The default will be set for maximum
##'     likelihood if \code{alpha==0.0}
##' @param nreps integer statnmg how many iterations to perform
##' @param ... unused
##' @return numeric vector with attributes \dQuote{"logLik"},
##'     \dQuote{"alpha"}, and \dQuote{"call"}
##' @export
##' @useDynLib ECTC
##' @importFrom Rcpp evalCpp
##' @author Charles Berry
##' 
##' @author Charles Berry
fixedEtaEM <-
    function(wtab,eta,params,alpha,p.init=NULL,
	     weights=NULL,nreps=1L,...)
{
    mc <- match.call()
    eodp <- eta %*% params$omega %*% diag(params$psi)
    eodp <- eodp/rowSums(eodp)
    lkMatrix <- t(dmulti(wtab$tab,eodp, .Machine$double.xmin))
    ritab <-  n <- wtab$n
    n.d <- ncol(lkMatrix)
    n.g <- nrow(lkMatrix)
    if (alpha==0.0){
	weights <- rep(1.0,nrow(lkMatrix))
	alpha.arg <- 1.0
    } else {
	if (is.null(weights)) weights <-
				  rep(1.0,nrow(lkMatrix))
	weights <- weights/sum(weights)
	alpha.arg <- alpha
    }
    prob.g <- if (!is.null(p.init)) p.init else weights
    stopifnot(length(prob.g)==nrow(lkMatrix))
    prob.g <- fixedEtaEMCall( lkMatrix,
                             prob.g, as.matrix(ritab), alpha.arg, weights,nreps ) 
    ## for (ir in 1:nreps){
    ##     prob.g.giv.d <- lkMatrix*prob.g
    ##     prob.g.giv.d <- prob.g.giv.d/rep(colSums(prob.g.giv.d),each=n.g)
    ##     k.d <- prob.g.giv.d %*% as.matrix(ritab)
    ##     prob.g <- prop.table(pmax(0,alpha*weights -1 + as.vector(k.d)))  
    ##    }
    dim(prob.g) <- NULL
    llk <- sum(n*log(prob.g %*% lkMatrix))
    ##    loglik <- dpLogLik(lkMatrix,prob.g,data.index,n)
    attr(prob.g,"logLik") <- llk
    attr(prob.g,"alpha") <- alpha
    attr(prob.g,"call") <- mc
    class(prob.g) <- "fixedEtaEM"
    prob.g}

fixedEtaEMCall <- function(lkMatrix, probg, ritab, alpha, weights, nreps){
    .Call('fixedEtaEMCall', lkMatrix, probg, ritab, alpha, weights, nreps)
}
