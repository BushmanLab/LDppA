##' Simulate the effect of sampling multinomials, subjecting them to
##' error, and random truncation.
##'
##' 
##' @title Simulate Cell Type Data
##' @param n list of vectors of replicate numbers
##' @param m list of cell numbers to sample
##' @param eta matrix whose rows are cell type probability vectors
##' @param omega confusion matrix for reassigning cell types
##' @param psi probability of cell type observation
##' @return matrix with \code{length(psi)} columns and
##'     \code{sum(unlist(n))} or fewer rows
##' @export
##' @examples
##' simW(1,10,matrix(1:3,nrow=1)/6,diag(3),rep(1,3))
##' simW(list(rep(1,10),2),
##' 	list( 10*(1:10),500),
##' 	prop.table(rbind(c(1,4,4),c(4,4,1)),1),
##' 	diag(3), rep(0.5,3))
##' @author Charles Berry
simW <- function(n,m,eta,omega,psi)
{
    stopifnot(length(n)==length(m),length(m)==nrow(eta))
    stopifnot(lengths(n)==lengths(m))
    stopifnot(all.equal(rowSums(eta),rep(1.0,nrow(eta))))
    n.vec <- unlist(n,use.names=FALSE)
    m.orig <- m.rep <- rep( unlist(m,use.names=FALSE), n.vec)
    m.len <- length(m.rep)
    eta.row <- rep(1:nrow(eta),sapply(as.list(n),sum))
    k <- length(psi)
    eop <- eta %*% omega %*% cbind(diag(psi), 1-psi)
    p.continue <- eop / ( 1- eop %*% upper.tri(diag(k+1)))
    ## play it safe here:
    p.continue[p.continue>1.0 | is.na(p.continue)] <- 1.0
    x <- array(0.0,c(length(m.rep),k))
    for (i in 1:k){
        r <- rbinom(m.len,m.rep,p.continue[eta.row,i])
        x[,i] <- r
        m.rep <- m.rep - r
    }

    x[m.rep!=m.orig,]
}
