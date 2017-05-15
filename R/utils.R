##' Discrete Representation of Posterior
##'
##' The posterior is discretized by assigning the mass to the nearest
##' point on a latttice on a k simplex.
##' @title Discrete Lattice Mass
##' @param eta M by \code{k+1} numeric matrix of composition
##'     probabilities
##' @param z integer or numeric vector of M
##'     probabilities or counts giving the amount of mass associated
##'     with each row of \code{eta}.
##' @param V numeric vector of M stick-breaking
##'     probabilities to be converted to a probability like \code{z}.
##' @param nintervals integer value giving number of intervals in each
##'     dimension.
##' @param digits integer giving number of trailing digits to use in
##'     labels.
##' @return vector of probabilities named to indicate the lattice points.
##' @export
##' @importFrom stats xtabs
##' @author Charles Berry
massGrid <- function(eta,z,V,nintervals=10,digits=NULL){
    if (is.null(digits))
        digits <- if (nintervals%in%c(10,5,2)) 1 else 2
    missing.z <- missing(z)
    missing.V <- missing(V)
    if (missing.z && missing.V) stop("Must specify either z or V argument")
    k <- ncol(eta)-1
    lattice.levels <-
        do.call(paste,as.data.frame(t(simplexGrid(0,nintervals,k)/nintervals)))
    factor.eta <-
        factor(do.call(paste,as.data.frame(roundSimplex(eta,1/nintervals))),
	       lattice.levels)
    pr <- if (missing.V) { # relative frequencies
	      if ( z[ length(z) ] == 1.0 && max(z)==1.0)
                  warning("Second arg assumed to be z, but last element is 1.0")
	      prop.table(z)
          } else {      # stick breaking probabilities
	      dZ.V(V)
          }

    tab <- xtabs(pr~factor.eta)
    fmt <- paste0(" %.",digits,"f")
    rowvals <- do.call(rbind,lapply(strsplit(names(tab)," "),
                                    function(x) sprintf(fmt,as.numeric(x))))
    names(tab) <- do.call(paste,as.data.frame(rowvals))
    unclass(tab)
}

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
