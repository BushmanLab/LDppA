##' Maximum Likelihood Estimation for the ECTC Model 
##'
##' Under the ECTC model given starting values for eta and V, expected
##' values for \code{sum(Z==t)} and cell type counts given \code{Z}
##' can be obtained.  The rows of eta then can be updated using
##' maximum likelihood.
##' @title estimateMaxLik
##' @param V numeric vector of initial values for \code{V}. Last value
##'     is 1.0.
##' @param eta numeric matrix of initial values for \code{eta}.
##' @param alpha parameter for a Dirichlet prior on \code{dZ.V(V)} 
##' @param params list of hyperparameters with elements \code{omega},
##'     and \code{psi}. Additional elments are ignored.
##' @param tab the result of \code{\link{wttab}()}
##' @param max.iter \code{integer} limiting iteration
##' @param rel.step \code{numeric} value limiting iteration
##' @param abs.step \code{numeric} value limiting iteration
##' @param alpha.max upper limit of \code{alpha/nrow(eta)}
##' @export
##' @importFrom stats optim optimize dmultinom
##' @return \code{list} with elements \code{logLik}, \code{eta},
##'     \code{prob.z}, \code{V}, \code{iter}, \code{dllk} and \code{call}
##' @author Charles Berry
estimateMaxLik <-
    function(
	     V,eta,alpha=0,params,tab,max.iter=500L,
	     rel.step=1e-06,abs.step=1e-3,alpha.max=1.0)
{
    mc <- match.call()
    argmax.llk <- function(phi,w,omega.psi,tol=1e-10){
	log.eta <- c(0,phi)-max(0,phi)
	eta <- exp(log.eta)/sum(exp(log.eta))
	eta.op <- eta%*%omega.psi
	eta.op <- eta.op/sum(eta.op)
	res <- sum(w*log(eta.op))
	res
    }
    dllk <- function(phi,w,omega.psi){
	xp2 <- c(0,phi)
	xp3 <- exp(xp2)
	xp4 <- xp3/sum(xp3)
	xp5 <- xp4%*%omega.psi
	xp6 <- xp5/sum(xp5)
	dxp2.dphi <- rbind(0,diag(nrow=length(phi)))
	dxp3.dxp2 <- diag(xp3)
	dxp4.dxp3 <- diag(nrow=length(xp3))/sum(xp3) - xp3/sum(xp3)^2
	dxp5.dxp4 <- t(omega.psi)
	dxp6.dxp5 <- diag(nrow=length(xp5))/sum(xp5) - as.vector(xp5)/sum(xp5)^2
	dxp7.dxp6 <- diag(1/as.vector(xp6))
	dxp8.dxp7 <- matrix(w,nrow=1)
	dxp8.dxp7 %*% dxp7.dxp6 %*% dxp6.dxp5 %*% dxp5.dxp4 %*%
	    dxp4.dxp3 %*% dxp3.dxp2 %*% dxp2.dphi}
    opt.fun <- function(i){
	## using good starting values helps speed and accuracy
	## using dumb starting values gives some negative updates
	safe.eta <- prop.table(eta[i,]+1e-08)
	opt <- optim(log(safe.eta[-1])-log(safe.eta[1]),
		     argmax.llk,dllk,
		     w=y[i,],
		     omega=omega%*%diag(psi),
		     control=list(fnscale=-1))
	eta.from.phi(opt$par)
    }
    update.prob.z.w <-
	function(prob.z,lik.zw)
	    prop.table(lik.zw * prob.z, 2) # equiv diag(prob.z) %*% lik.zw
    update.prob.z <-
	function(prob.z.w,a=alpha,kv=k,wt=tab)
    {
	ez.w <- as.vector( prob.z.w %*% wt$n + a/kv - 1)
	prop.table(pmax(0,ez.w))
    }
    ldmn <- function(alpha,x){
	k <- length(x)
	n <- sum(x)
	lgamma(k*alpha) + lfactorial(n) - lgamma(n+k*alpha) +
	    sum(lgamma(x+alpha)) - sum(lfactorial(x)) - k*lgamma(alpha)
    }
    update.alpha <-
	function(ex){
	    optimize(
		function(x) ldmn(x,ex),
		c(0.0,alpha.max),maximum=TRUE)
	}
    omega <- params$omega
    psi <- params$psi
    k <- nrow(eta)
    ## inits
    eodp <- eta %*% omega %*% diag(psi)
    eodcp <- eta %*% omega %*% diag(1-psi)
    eop <- rowSums(eodp)
    eodp <- eodp/eop
    like.zw <- t(dmulti(tab$tab,eodp))
    prob.z <- dZ.V(V)
    if (alpha!=0) {
	prob.z.w <- update.prob.z.w(prob.z,like.zw)
	prob.z <- update.prob.z(prob.z.w)
	ex <- prob.z.w %*% tab$n
	alpha <- update.alpha(prob.z)$maximum * k # like Dir(alpha/N,...)
    }
    ## updates
    dllk <- Inf
    llk <- -Inf
    iter <- 0L
    while (iter<max.iter && dllk >= min(abs.step,abs(rel.step * llk)))
    {
	iter <- iter + 1L
	old.llk <- llk
	prob.z.w <- update.prob.z.w(prob.z,like.zw)
	prob.z <- update.prob.z(prob.z.w)
	if (alpha!=0) {
	    ex <- prob.z.w %*% tab$n
	    uap <- update.alpha(ex)
	    alpha <- uap$maximum * k # like Dir(alpha/N,...)
	    alphallk <- uap$objective
	}
	y <-  prob.z.w %*% (tab$tab*tab$n) 
	rowvals <- sapply(1:nrow(y), opt.fun)
	eta <- t(rowvals)
	eodp <- eta %*% omega %*% diag(psi)
	eodcp <- eta %*% omega %*% diag(1-psi)
	eop <- rowSums(eodp)
	eodp <- eodp/eop
	like.zw <- t(dmulti(tab$tab,eodp))
	llk <- sum(tab$n*log(prob.z%*%like.zw))
	dllk <- llk-old.llk
    }
    if (alpha==0) alphallk <- 0
    V <- prob.z/rev(cumsum(rev(prob.z)))
    list(logLik=llk,eta=eta,alpha=alpha,prob.z=prob.z,V=V,iter=iter,
	 dllk=llk-old.llk,alpha.llk=alphallk,call=mc)
}

##' MCMC Sampler for ECTC Model
##'
##' Under the ECTC model with subsampling given eta and V, z can be
##' drawn collapsing over X and Y, then full Conditionals for V,
##' alpha, eta, and X are drawn. The draws for X are are from
##' multinomial distributions and many can be combined as they share
##' the same probability vector and will be summed to update eta. The
##' draws for Y can collapse over X and R.
##' @title estimateComps
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
##' @param fix.eta if \code{TRUE} do not update eta
##' @param prob.Z.dirich if \code{TRUE} sample composition probability
##'     by Dirichlet using \code{rep(alpha,nrow(eta))} as the
##'     parameter
##' @return \code{list} with components \code{monitors} (recording
##'     components of the log posterior), \code{Vs} (values of V),
##'     \code{etas} (values of eta), and \code{call} (the call).
##' @export
##' @importFrom stats rbeta rgamma rmultinom dnbinom rnbinom
##' @importFrom utils head
##' @useDynLib ECTC
##' @author Charles Berry
estimateComps <- function(V,eta,alpha,params,tab, nreps=1L, nburn=0L, nthin=1L,
                          save.last=FALSE,fix.eta=FALSE,prob.Z.dirich=FALSE)
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
        loop <- .C("sampleCTC",
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
                   xsums=integer(T*ka),
                   fixeta=as.integer(fix.eta),
                   pzdirich=as.integer(prob.Z.dirich))
        eta <- matrix(loop$eta,nrow=T)
        V <- loop$V
        alpha <- loop$alpha
    }


    for (i in 1:nreps){

        loop <- .C("sampleCTC",
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
                   xsums=integer(T*ka),
                   fixeta=as.integer(fix.eta),
                   pzdirich=as.integer(prob.Z.dirich))

        zy.tab <- matrix(loop$zy,nrow=T)
        eta <- matrix(loop$eta,nrow=T)
        V <- loop$V
        alpha <- loop$alpha
        zsums <- rowSums(zy.tab)
        ## monitors
        vals <- c(
            alpha= logP.alpha(alpha,s),
            V= logP.V.alpha(V,alpha),
            z= logP.ztab.v(zsums,V),
            eta = logP.eta(eta,lamb), ## (lambda == 1) ==> constant
            Y = logP.Y.eta.V(W,eta,V,omega,n,psi,eps)
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
