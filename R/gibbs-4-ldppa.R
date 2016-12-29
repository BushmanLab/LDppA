##' MCMC Sampler for LDppA Model - fixed eta
##'
##' Under the LDppA model with subsampling given eta and V, z can be
##' drawn collapsing over X and Y, then full Conditionals for V,
##' alpha, and etaare drawn. Here no draws for X or eta are done. That
##' is, the allowable values of eta are pre-fixed.
##' @title ldppa.gibbs.4
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
##' @param nreps number of iteration groups to run after burn-in
##' @param nburn number of iterations to run before retaining samples
##' @param nthin iterate this many times between retaining values
##' @param save.last logical.  If TRUE save all of the last iteration.
##' @return \code{list} with components \code{monitors} (recording
##'     components of the log posterior), \code{Vs} (values of V),
##'     \code{etas} (values of eta), and \code{call} (the call).
##' @export
##' @importFrom stats rbeta rgamma rmultinom dnbinom rnbinom
##' @importFrom utils head
##' @useDynLib LDppA
##' @author Charles Berry
ldppa.gibbs.4 <- function(V,eta,alpha,params,tab, nreps=1L, nburn=0L, nthin=1L,
                          save.last=FALSE)
{
    sc <- match.call()
    ## params
    omega <- params$omega
    lamb <- params$lambda
    s <- params$s
    eps <- params$epsilon
    psi <- params$psi
    K <- ncol(tab$tab)
    ka <- ko <- K
    T <- nrow(eta)
    ## `y' is used below where `W' might have better mnemonic
    ## meaning. This is a legacy of code written before subsampling was
    ## taken into account.
    W <- tab$tab
    Wplus <- rowSums(W)
    Nw <- length(Wplus)
    n <- tab$n
    ndat <- as.integer(nrow(W))
    ## calc once
    omPsi <- omega%*%psi
    omDPsi <- omega%*%diag(psi)
    omCPsi <- omega%*%(1-psi)

    ## initialize
    monitors <- list()
    etas <- list()
    Vs <- list()
    alphas <- list()
    zs <- list()

    ikeep <- 0

    if (nburn>0){
        loop <- .C("fixedEta",
                   reps=as.integer(nburn),
                   T=as.integer(T),
                   ka=as.integer(K),
                   ko=as.integer(K),
                   n=as.integer(n),
                   ndat=as.integer(ndat),
                   w=as.integer(W),
                   wp=as.integer(Wplus),
                   s=as.double(s),
                   lamb=as.double(lamb),
                   omcp=as.double(omCPsi),
                   omdp=as.double(omDPsi),
                   eps=as.double(eps),
                   V=as.double(V),
                   alpha=as.double(alpha),
                   eta=as.double(eta),
                   zy=integer(T*ndat),
                   eoy=double(ka),
                   etaomdp=double(T*ka),
                   prw=double(T),
                   pz=double(T),
                   workT=double(T),
                   xstmp=integer(ka),
                   xsums=integer(T*ka))
        eta <- matrix(loop$eta,nrow=T)
        V <- loop$V
        alpha <- loop$alpha
    }


    for (i in 1:nreps){

        loop <- .C("fixedEta",
                   reps=as.integer(nthin),
                   T=as.integer(T),
                   ka=as.integer(K),
                   ko=as.integer(K),
                   n=as.integer(n),
                   ndat=as.integer(ndat),
                   w=as.integer(W),
                   wp=as.integer(Wplus),
                   s=as.double(s),
                   lamb=as.double(lamb),
                   omcp=as.double(omCPsi),
                   omdp=as.double(omDPsi),
                   eps=as.double(eps),
                   V=as.double(V),
                   alpha=as.double(alpha),
                   eta=as.double(eta),
                   zy=integer(T*ndat),
                   eoy=double(ka),
                   etaomdp=double(T*ka),
                   prw=double(T),
                   pz=double(T),
                   workT=double(T),
                   xstmp=integer(ka),
                   xsums=integer(T*ka))

        zy.tab <- matrix(loop$zy,nrow=T)

        x.sums <- matrix(loop$xsums,nrow=T)

        eta <- matrix(loop$eta,nrow=T)

        V <- loop$V

        alpha <- loop$alpha

        zsums <- rowSums(zy.tab)

        ## monitor
        vals <- c(
            alpha= logP.alpha(alpha,s),
            V= logP.V.alpha(V,alpha),
            z= logP.ztab.v(zsums,V),
            eta = logP.eta(eta,lamb), ## (lambda == 1) ==> constant
            Y = logP.Y.eta.V(W,eta,V,omega,n,psi)
        )


        monitors[[i]] <- vals
        etas[[i]] <- eta
        Vs[[i]] <- V
        alphas[[i]] <- alpha
        zs[[i]] <- zsums
    }


    list(monitors=monitors,etas=etas,Vs=Vs,alphas=alphas,zs=zs,call=sc,
         last=if (save.last) loop else NULL )
}
