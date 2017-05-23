##' @importFrom stats dbeta dgamma

logP.eta <- function(eta,lambda){
    sum(ddirichlet(eta,lambda,log.p=TRUE))
}

log.prob.v <- function(v){
    v <- pmin(1,v)
    log(c(v,1))+c(0,cumsum(log1p(-v)))
}

logP.V.alpha <- function(V,alpha){
    sum(dbeta(head(V,-1),1,alpha,log=TRUE))
}

logP.alpha <- function(alpha,s){
    dgamma(alpha,s[1],rate=s[2],log=TRUE)
}

logP.Y.eta.V <- function(y,eta,V,omega,n=1,psi=NULL,epsilon=c(0,0)){
    prob <- exp(log.prob.v(V[-length(V)]))
    pr.y.eta.t <- eta%*%omega
    pr.seen <- rowSums(pr.y.eta.t)
    if (!is.null(psi)){
        pr.y.eta.t <- sweep(pr.y.eta.t,2,psi,"*")
    }
    pr.multi <- dmulti(y,quick.prop.table(pr.y.eta.t,1))
    if (sum(epsilon)>0){
        pr.wp <-
            outer(rowSums(y),pr.seen,
                  function(x,y) dnbinom(x,epsilon[1],
                                        1-pr.seen/(pr.seen+epsilon[2])))
        pr.multi <- pr.multi*pr.wp
    }
    llk <- pr.multi %*% prob
    sum(n*log(llk))
}
