##' Mixture Component Probability from Stick Breaking
##'
##' Given a collection of values, \code{V}, of fractions (except the
##' last which is 1), from a stick breaking process, find the
##' probability of each class.
##' @title Stick Breaking Probabilities
##' @param V vector of values in the unit interval except
##'     \code{V[length(V)]==1}
##' @return a vector of probabilities or a matrix whose rows contain
##'     such vectors
##' @export
##' @author Charles Berry
dZ.V <- function(V){
    if (is.list(V)) V <- do.call(rbind,V)
    if (is.matrix(V)){
        nc <- ncol(V)
        V* exp( log1p(-V[,-nc])%*%upper.tri(diag(nc))[-nc,] )
    } else {
        prob.z.v(V)
    }
}

prob.z.v <- function(v) v*exp(c(0,cumsum(log1p(-v)[-length(v)])))

##' @importFrom stats dmultinom
logP.ztab.v <- function(ztab,v){
    pie <- prob.z.v(v)
    ## unvectorized dmultinom() is OK here
    dmultinom(ztab,prob=pie,log=TRUE)
}
