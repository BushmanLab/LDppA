##' @importFrom stats cutree hclust dist coef lm optim
gamma.initial <- function(Ntau) cbind(1,(Ntau-1):1)

## this can be hastened by replacing prop.table() 
mu.initial <- function(tab,omega){
    tab <- tab$tab
    K <- ncol(omega)
    N <- nrow(tab)
    res <-
        array(rep(tab,K)*rep(omega,each=N), c(N,K,K))
    res <- quick.prop.table(res+1/K^2,c(1,3))
    res}

tau.initial <-
    function(tab,omega,lambda,k=NULL,max.rows=NULL) {
        eta.from.phi <-
            function(phi) c(1,exp(phi))/(1+sum(exp(phi)))
        argmax.llk <- function(phi,w,omega.psi,tol=1e-10){
            eta <- c(1,exp(phi))/(1+sum(exp(phi)))
            eta.op <- eta%*%omega.psi
            eta.op <- eta.op/sum(eta.op)
            res <- dmultinom(w,prob=eta.op,log=TRUE)
            res
        }
        y <- tab$tab
        wt <- tab$n
        cells <- rowSums(y)*wt
        if (is.null(max.rows)) max.rows <-
                                   min(sum(wt),k*100)
        if (is.null(k)) k <- length(lambda)
        boot.samp <- sample(nrow(y),max.rows,replace=TRUE,prob=wt)
        rowvals <- apply(y,1,
                         function(x){
                             opt <- optim(rep(0,length(lambda)-1),
                                          argmax.llk,
                                          w=x,
                                          omega=omega,
                                          control=list(fnscale=-1))
                             eta.from.phi(opt$par)
                         })

        y <- t(rowvals[,boot.samp])*cells[boot.samp]
        ind <- cutree(hclust(dist(y)),k=k)
        unname(coef(lm(y~0+as.factor(ind)))+rep(lambda,each=k))
    }

x.n.plus.r.plus <- function(y,mu)
    rowSums(mu*
            as.vector(y[,rep(1:ncol(y),each=ncol(y))]),dims=2)

Snt.fun <- function(y,mu,gamma,tau){
    xnr <- x.n.plus.r.plus(y,mu)
    ElogV <- c(digamma(gamma[,1])-digamma(rowSums(gamma)),0)
    sumElogcV <- cumsum(c(0,digamma(gamma[,2])-digamma(rowSums(gamma))))
    dtau <- digamma(tau)-digamma(rowSums(tau))
    rep(ElogV+sumElogcV,each=nrow(mu))+
        sapply(1:nrow(tau),function(tee) xnr%*%dtau[tee,])
}

rowMax <- function(x) do.call(pmax,as.data.frame(x))

phi.initial <- function(y,mu,gamma,tau){
    Snt <- Snt.fun(y,mu,gamma,tau)
    quick.prop.table(exp(Snt-rowMax(Snt)),1)
}
