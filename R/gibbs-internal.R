##' @importFrom stats dbeta dgamma

logP.eta <- function(eta,lambda){
    sum(ddirichlet(eta,lambda,log.p=TRUE))
}

log.prob.v <- function(v) log(c(v,1))+c(0,cumsum(log1p(-v)))

logP.V.alpha <- function(V,alpha){
    sum(dbeta(head(V,-1),1,alpha,log=TRUE))
}

logP.alpha <- function(alpha,s){
    dgamma(alpha,s[1],rate=s[2],log=TRUE)
}

logP.Y.eta.V <- function(y,eta,V,omega,n=1){
    prob <- exp(log.prob.v(V[-length(V)]))
    pr.y.eta.t <- eta%*%omega
    llk <- dmulti(y,pr.y.eta.t)%*%prob
    sum(n*log(llk))
}
