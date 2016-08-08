##' MCMC Sampler for LDppA Model
##'
##' Under the LDppA model given eta and V, z can be drawn collapsing
##' over X, then full Conditionals for V, alpha, eta, and X are
##' drawn. The draws for X are are from multinomial distributions and
##' many can be combined as they share the same probability vector and
##' will be summed to update eta.
##' @title ldppa.gibbs
##' @param V numeric vector of initial values for \code{V}.
##' @param eta numeric matrix of initial values for \code{eta}.
##' @param alpha numeric value to initialize \code{alpha}.
##' @param params list of hyperparameters with elements \code{omega},
##'     \code{lambda}, and \code{s}
##' @param tab the result of \code{\link{wttab}()}
##' @param nreps number of iterations to run after burn-in
##' @param nburn number of iterations to run before retaining samples
##' @param nthin retain values for \code{nreps/nthin} of the
##'     post-burn-in iterations
##' @return \code{list} with components \code{monitors} (recording
##'     components of the log posterior), \code{Vs} (values of V},
##'     \code{etas} (values of eta), and \code{call} (the call).
##' @export
##' @importFrom stats rbeta rgamma rmultinom
##' @importFrom utils head
##' @author Charles Berry
ldppa.gibbs <- function(V,eta,alpha,params,tab, nreps=1L, nburn=0L, nthin=1L)
{
    sc <- match.call()
    ## params
    omeg <- params$omega
    lamb <- params$lambda
    s <- params$s
    K <- ncol(tab$tab)
    T <- nrow(eta)
    y <- tab$tab
    n <- tab$n
    ## initialize
    monitors <- list()
    etas <- list()
    Vs <- list()

    ikeep <- 0
    for (i in 1:(nreps+nburn)){
        ## sample z
        logprob.kernel <- tcrossprod(log(eta%*%omeg),y)
        prob.z.giv.y.v.eta <-
    	prop.table(
                exp(sweep(logprob.kernel,2,
                          apply(logprob.kernel,2,max),"-")) *
                prob.z.v(V),2)
        zy.tab <- sapply(1:length(n),
                         function(x) rmultinom(1,n[x],prob.z.giv.y.v.eta[,x]))
        ## z.list <-
     ##     lapply(1:length(n),
        ##            function(x) sample.int(T,n[x],
        ##                                   prob=prob.z.giv.y.v.eta[,x],
        ##                                   replace=TRUE))

        ## sample x

        eo.y <-
            sapply(1:K,function(x) {
                eo <- quick.prop.table(sweep(eta,2,omeg[,x],"*"),1)
                eo})
        dim(eo.y) <- c(T,K,K)

        y.z <- zy.tab %*% y

        x.sums <- matrix(0.0,nrow=T,ncol=K)

        for (t in 1:T)
            for (k in 1:K)
                x.sums[t,] <- x.sums[t,] + rmultinom(1,y.z[t,k],eo.y[t,,k])

        ## sample eta
        lamb.new <- sweep(x.sums,2,lamb,"+")
        eta <- rdirichlet(lamb.new)
        ## sample V
        z.tab <- rowSums(zy.tab)
        z.gt <- sum(z.tab)-cumsum(z.tab)
        V <- c(rbeta( length(V)-1, 1+head(z.tab,-1),
    		 alpha+head(z.gt,-1)),
    	   1)
        ## sample alpha
        alpha <- rgamma(1,s[1]+T-1,rate=s[2]-sum(log(1-head(V,-1))))
        ## monitor
        vals <- c(
    	alpha= logP.alpha(alpha,s),
    	V= logP.V.alpha(V,alpha),
    	z= logP.ztab.v(z.tab,V),
    	eta = logP.eta(eta,lamb), ## (lambda == 1) ==> constant
    	Y = logP.Y.eta.V(y,eta,V,omeg,n)
        )

        if ( i>nburn &&
             ((i-nburn-1)%%nthin==0 || i == (nreps+nburn) ))
        {
    	ikeep <- ikeep +1
    	monitors[[ikeep]] <- vals
    	etas[[ikeep]] <- eta
    	Vs[[ikeep]] <- V
        }
    }

    list(monitors=monitors,etas=etas,Vs=Vs,call=sc)
}
