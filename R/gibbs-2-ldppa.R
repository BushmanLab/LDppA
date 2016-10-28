##' MCMC Sampler for LDppA Model
##'
##' Under the LDppA model with subsampling given eta and V, z can be
##' drawn collapsing over X and Y, then full Conditionals for V,
##' alpha, eta, and X are drawn. The draws for X are are from
##' multinomial distributions and many can be combined as they share
##' the same probability vector and will be summed to update eta. The
##' draws for Y can collapse over X and R.
##' @title ldppa.gibbs.2
##' @param V numeric vector of initial values for \code{V}. Last value
##'     is 1.0.
##' @param eta numeric matrix of initial values for \code{eta}.
##' @param alpha numeric value to initialize \code{alpha}.
##' @param params list of hyperparameters with elements \code{omega},
##'     \code{lambda}, \code{psi}, \code{epsilon}, and \code{s}. Note
##'     that \code{epsilon[1]} must be an integer, but can be set to
##'     zero (which with \code{epsilon[2]==0} corresponds to a
##'     non-informative prior for the number of cells sampled).
##' @param tab the result of \code{\link{wttab}()}
##' @param nreps number of iterations to run after burn-in
##' @param nburn number of iterations to run before retaining samples
##' @param nthin retain values for \code{nreps/nthin} of the
##'     post-burn-in iterations
##' @return \code{list} with components \code{monitors} (recording
##'     components of the log posterior), \code{Vs} (values of V),
##'     \code{etas} (values of eta), and \code{call} (the call).
##' @export
##' @importFrom stats rbeta rgamma rmultinom dnbinom rnbinom
##' @importFrom utils head
##' @author Charles Berry
ldppa.gibbs.2 <- function(V,eta,alpha,params,tab, nreps=1L, nburn=0L, nthin=1L)
{
    sc <- match.call()
    ## params
    omega <- params$omega
    lamb <- params$lambda
    s <- params$s
    eps <- params$epsilon
    psi <- params$psi
    K <- ncol(tab$tab)
    T <- nrow(eta)
    ## `y' is used below where `W' might have better mnemonic
    ## meaning. This is a legacy of code written before subsampling was
    ## taken into account.
    W <- tab$tab
    Wplus <- rowSums(W)
    Nw <- length(Wplus)
    n <- tab$n
    ## calc once
    omPsi <- omega%*%psi
    omDPsi <- omega%*%diag(psi)
    omCPsi <- omega%*%(1-psi)

    ## initialize
    monitors <- list()
    etas <- list()
    Vs <- list()
    alphas <- list()

    ikeep <- 0
    for (i in 1:(nreps+nburn)){
        ## z
        ## integrate out R
        ## W|W[+],Z is multinomial
        ## W[+]|Z is neg-binom

        etaOmDPsi <- eta %*% omDPsi
        ## prob observe any event given t
        PrW <- rowSums(etaOmDPsi) 

        etaOmCPsi <- quick.prop.table(sweep(eta,2,omCPsi,"*"),1)

        ## sample z
        ## all conditioned on Z:
        logPrWplus <-
            if (eps[1]==0)
                0
            else 
                outer(PrW,Wplus, function(p,N)
                    dnbinom( N, eps[1], 1-eps[2]/p, log=TRUE))
        logPrW.giv.Wplus <-
            tcrossprod(log(etaOmDPsi/PrW),W)
        logPrW <- logPrWplus+logPrW.giv.Wplus

        pr.Z.giv.W.V.eta <-
            prop.table(
                exp(sweep(logPrW,2,apply(logPrW,2,max),"-")) *
                prob.z.v(V) ,2)

        zy.tab <- sapply(1:length(n),
                         function(x) rmultinom(1,n[x],pr.Z.giv.W.V.eta[,x]))
        ## sample x
        ## integrate out R
        ## X[-]|W[+],Z is neg-binom
        ## X|W,Z is multinomial

        eo.y <-
            sapply(1:K,function(x) {
                eo <- quick.prop.table(sweep(eta,2,omDPsi[,x],"*"),1)
                eo})
        dim(eo.y) <- c(T,K,K)

        y.z <- zy.tab %*% W

        x.sums <- matrix(0.0,nrow=T,ncol=K)

        for (t in 1:T){
            for (k in 1:K)
                x.sums[t,] <- x.sums[t,] + rmultinom(1,y.z[t,k],eo.y[t,,k])
            zy.pos <- zy.tab[t,]>0
            ## (1-prob seen)/(1+rate parm for prior)
            nb.parm <- (1-PrW[t])/(1+eps[2])

            x.unseen <-
                rnbinom(sum(zy.pos),
                        zy.tab[t,zy.pos] * (eps[1]+Wplus[zy.pos]),
                        1-nb.parm)

            x.sums[t,] <- x.sums[t,]+
                rmultinom(1,sum(x.unseen),etaOmCPsi[t,])
        }

        ## sample eta
        lamb.new <- sweep(x.sums,2,lamb,"+")
        eta <- rdirichlet(lamb.new)
        ## sample V
        z.tab <- rowSums(zy.tab)
        z.gt <- sum(z.tab)-cumsum(z.tab)
        ## ensure V[-T] < 1.0
        V <- c(
            pmin(0.99, rbeta( length(V)-1, 1+head(z.tab,-1),
                             alpha+head(z.gt,-1))),
            1)
        ## sample alpha
        alpha <- rgamma(1,s[1]+T-1,rate=s[2]-sum(log(1-head(V,-1))))

        ## monitor
        vals <- c(
            alpha= logP.alpha(alpha,s),
            V= logP.V.alpha(V,alpha),
            z= logP.ztab.v(z.tab,V),
            eta = logP.eta(eta,lamb), ## (lambda == 1) ==> constant
            Y = logP.Y.eta.V(W,eta,V,omega,n,psi)
        )

        if ( i>nburn &&
             ((i-nburn-1)%%nthin==0 || i == (nreps+nburn) ))
        {
            ikeep <- ikeep +1
            monitors[[ikeep]] <- vals
            etas[[ikeep]] <- eta
            Vs[[ikeep]] <- V
            alphas[[ikeep]] <- alpha
        }
    }

    list(monitors=monitors,etas=etas,Vs=Vs,alphas=alphas,call=sc)
}
