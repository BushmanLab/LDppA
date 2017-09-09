##' Maximum Likelihood or Posterior Estimation for the ECTC Model 
##'
##' Under the ECTC model given starting values for eta and V, expected
##' values for \code{sum(Z==t)}.  Maximum likelihood or maximum
##' posterior estimates of cell type counts given \code{Z} and of eta
##' then can be obtained.  The value of \code{alpha} (if not
##' initialized to zero) is obtained by maximizing the posterior.  As
##' such, this qualifies as an empirical Bayes approach.
##' @title estimateMaxLik
##' @aliases estimateMaxLik2
##' @usage estimateMaxLik(V,eta,alpha=0,params,tab,max.iter=500L,
##'        rel.step=1e-06,abs.step=1e-5,alpha.max=1.0,...)
##'
##' estimateMaxLik2(V,eta,alpha=0,params,tab,max.iter=500L,
##'        rel.step=1e-06,abs.step=1e-5,alpha.max=1.0,...)
##' @param V numeric vector of initial values for \code{V}. Last value
##'     is 1.0.
##' @param eta numeric matrix of initial values for \code{eta}.
##' @param alpha If \code{alpha==0.0} maximum likelihood estimates
##'     will be obtained. if non-zero, \code{rep(alpha/nrow(eta),
##'     nrow(eta))} is the initial value of the parameter for a
##'     Dirichlet prior on \code{prob.z} (aka \code{dZ.V(V)}).
##' @param params list of hyperparameters with elements \code{omega},
##'     and \code{psi}. Additional elments are ignored.
##' @param tab the result of \code{\link{wttab}()}
##' @param max.iter \code{integer} limiting iteration
##' @param rel.step \code{numeric} value limiting iteration
##' @param abs.step \code{numeric} value limiting iteration
##' @param alpha.max upper limit of \code{alpha/nrow(eta)}
##' @param ... unused
##' @export
##' @importFrom stats optim optimize dmultinom
##' @return \code{list} with elements \code{logLik}, \code{eta},
##'     \code{prob.z}, \code{V}, \code{iter}, \code{dllk} and
##'     \code{call}
##' @author Charles Berry
estimateMaxLik <-
    function(
	     V,eta,alpha=0,params,tab,max.iter=500L,
	     rel.step=1e-06,abs.step=1e-5,alpha.max=1.0,...)
{
    mc <- match.call()
    argmax.llk <- function(phi,w,omega.psi){
	.Call("amllk",phi,as.double(w),omega.psi)
    }
    dllkdphi <- function(phi,w,omega.psi)
	.Call("dldphi",
	      as.double(phi), as.double(w), as.double(omega.psi))
    opt.fun <- function(i){
	## using good starting values helps speed and accuracy
	## using dumb starting values gives some negative updates
	safe.eta <- prop.table(eta[i,]+1e-08)
	opt <- optim(log(safe.eta[-1])-log(safe.eta[1]),
		     argmax.llk,dllkdphi,
		     w=y[i,],
		     omega=omega%*%diag(psi),
		     method="BFGS",
		     control=list(fnscale=-1))
	eta.from.phi(opt$par)
    }
    all.derivs <- function(y,eta,omega,psi){
	op <- omega%*%diag(psi)
	safe.eta <- prop.table(eta+1e-08,1)
	log.safe.eta <- log(safe.eta[,-1,drop=FALSE])-log(safe.eta[,1])
	ne <- nrow(eta)
	ny <- nrow(y)
	sapply(1:ne, function(i){
	    dde <- sapply(1:ny,
			  function(j) dllkdphi(log.safe.eta[i,],y[j,],op))
	    rowSums(dde)})}
    update.prob.z.w <-
	function(prob.z,lik.zw){
	    x <- lik.zw * prob.z
	    x/rep(colSums(x),each=length(prob.z))
        }
    update.prob.z <-
	function(prob.z.w,a=alpha,kv=k,wt=tab)
    {
	ez.w <- as.vector( prob.z.w %*% wt$n + a/kv - if (a==0) 0 else 1 )
	x <- pmax(1e-12,ez.w)
	x/sum(x)
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
    like.zw <- t(dmulti(tab$tab,eodp, .Machine$double.xmin))
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
    while (iter<max.iter && (
	dllk >= min(abs.step,abs(rel.step * llk)) ||
	max(abs(old.prob.z-prob.z)>abs.step)))
    {
	iter <- iter + 1L
	old.llk <- llk
	prob.z.w <- update.prob.z.w(prob.z,like.zw)
	old.prob.z <- prob.z
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
	## there can be error here for low probabilities
	## as a check: completellk >= like.zw
	like.zw <- t(dmulti(tab$tab,eodp,.Machine$double.xmin))
	llk <- sum(tab$n*log(prob.z%*%like.zw))
	dllk <- llk-old.llk
    }
    ## complete data loglike
    ex.z.w <- sweep(prob.z.w,2,tab$n,"*")
    completellk <- sum(ex.z.w*t(dmulti(tab$tab, eodp, log.p=TRUE)))
    if (alpha==0) alphallk <- 0
    V <- prob.z/rev(cumsum(rev(prob.z)))
    list(logLik=llk,eta=eta,alpha=alpha,prob.z=prob.z,V=V,iter=iter,
	 dllk=llk-old.llk,alpha.llk=alphallk,complete.llk=completellk,
	 call=mc)
}

##' @export
estimateMaxLik2 <-
    function(
	     V,eta,alpha=0,params,tab,max.iter=500L,
	     rel.step=1e-06,abs.step=1e-5,alpha.max=1.0,...)
{
    mc <- match.call()
    argmax.llk <- function(phi,w,omega.psi){
	.Call("amllk",phi,as.double(w),omega.psi)
    }
    dllkdphi <- function(phi,w,omega.psi)
	.Call("dldphi",
	      as.double(phi), as.double(w), as.double(omega.psi))
    opt.fun <- function(i){
	## using good starting values helps speed and accuracy
	## using dumb starting values gives some negative updates
	safe.eta <- prop.table(eta[i,]+1e-08)
	opt <- optim(log(safe.eta[-1])-log(safe.eta[1]),
		     argmax.llk,dllkdphi,
		     w=y[i,],
		     omega=omega%*%diag(psi),
		     method="BFGS",
		     control=list(fnscale=-1))
	eta.from.phi(opt$par)
    }
    all.derivs <- function(y,eta,omega,psi){
	op <- omega%*%diag(psi)
	safe.eta <- prop.table(eta+1e-08,1)
	log.safe.eta <- log(safe.eta[,-1,drop=FALSE])-log(safe.eta[,1])
	ne <- nrow(eta)
	ny <- nrow(y)
	sapply(1:ne, function(i){
	    dde <- sapply(1:ny,
			  function(j) dllkdphi(log.safe.eta[i,],y[j,],op))
	    rowSums(dde)})}
    update.prob.z.w <-
	function(prob.z,lik.zw,prob.wp.z,wpui){
	    x <- (lik.zw * prob.wp.z[,wpui]) * prob.z 
	    x/rep(colSums(x),each=length(prob.z))
        }
    update.prob.z <-
	function(prob.z.w,a=alpha,kv=k,wt=tab)
    {
	ez.w <- as.vector( prob.z.w %*% wt$n + a/kv - if (a==0) 0 else 1 )
	x <- pmax(1e-12,ez.w)
	x/sum(x)
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
    wplus <- rowSums(tab$tab)
    wpu <- sort(unique(wplus))
    wpui <- match(wplus,wpu)
    ## mm <- model.matrix(~0+factor(wplus,wpu))
    mm <- diag(length(wpu))[wpui,]
    eodp <- eta %*% omega %*% diag(psi)
    eodcp <- eta %*% omega %*% diag(1-psi)
    eop <- rowSums(eodp)
    eodp <- eodp/eop
    like.zw <- t(dmulti(tab$tab,eodp, .Machine$double.xmin))
    prob.wp.z <- matrix(1.0,nrow=nrow(like.zw),ncol=length(wpu))
    prob.z <- dZ.V(V)
    if (alpha!=0) {
	prob.z.w <- update.prob.z.w(prob.z,like.zw,prob.wp.z,wpui)
	prob.z <- update.prob.z(prob.z.w)
	ex <- prob.z.w %*% tab$n
	alpha <- update.alpha(prob.z)$maximum * k # like Dir(alpha/N,...)
    }
    ## updates
    dllk <- Inf
    llk <- -Inf
    iter <- 0L
    while (iter<max.iter && (
	dllk >= min(abs.step,abs(rel.step * llk)) ||
	max(abs(old.prob.z-prob.z)>abs.step)))
    {
	iter <- iter + 1L
	old.llk <- llk
	prob.z.w <- update.prob.z.w(prob.z,like.zw,prob.wp.z,wpui)
	prob.wp.z <- prop.table(sweep(prob.z.w,2,tab$n,"*")%*%mm,1)
	old.prob.z <- prob.z
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
	## there can be error here for low probabilities
	## as a check: completellk >= like.zw
	like.zw <- t(dmulti(tab$tab,eodp,.Machine$double.xmin))
	llk <- sum(tab$n*log(prob.z%*%like.zw))
	dllk <- llk-old.llk
    }
    ## complete data loglike
    ex.z.w <- sweep(prob.z.w,2,tab$n,"*")
    completellk <- sum(ex.z.w*t(dmulti(tab$tab, eodp, log.p=TRUE)))
    if (alpha==0) alphallk <- 0
    V <- prob.z/rev(cumsum(rev(prob.z)))
    list(logLik=llk,eta=eta,alpha=alpha,prob.z=prob.z,V=V,iter=iter,
	 dllk=llk-old.llk,alpha.llk=alphallk,complete.llk=completellk,
	 call=mc)
}

##' Find Parsimonious MaxLik by Splitting and Merging
##'
##' This is a wrapper for \code{\link{estimateMaxLik}}, which see for
##' details about parameters and the return value
##' @title Stepwise Maximization
##' @param V sticks
##' @param eta compositions
##' @param alpha tuning
##' @param params omega and psi
##' @param tab \code{wttab(...)}
##' @param tol how fine to tune
##' @param max.cycles limit on number of cycles of merging/splitting 
##' @param ... other args to pass down
##' @param min.n smallest pseudo n for BIC
##' @param verbose print extra stuff, if 2 print lots
##' @return list with updated parms
##' @export
##' @author Charles Berry
stepMaxLik <-  function(V,eta,alpha,params,tab,tol=1e-6,max.cycles=10L,...,
			min.n=2,verbose=FALSE){
    try.again <- 2
    while (try.again && {max.cycles <- max.cycles-1} >= 0 ){
	mres <- mergeMaxLik(V,eta,alpha,params,tab,
			    min.n=min.n, ..., verbose=verbose)
	if (verbose==2) {
	    print(round(mres$eta,3))
	    print(round(mres$prob.z,4))
	}
	sres <- splitMaxLik(mres$V,mres$eta,mres$alpha,
			    params,tab,min.n=min.n,
			    ...,verbose=verbose)
	if (verbose==2) {
	    print(round(sres$eta,3))
	    print(round(sres$prob.z,4))
	}
	if (try.again==2)
	{
	    try.again <- 1
	    old.llk <- sres$logLik
	} else {
	    llk <- sres$logLik
	    try.again <- llk > old.llk+tol
	    old.llk <- llk
	}
	if (verbose) cat(old.llk,"\n")
	V <- sres$V
	eta <- sres$eta
	alpha <- sres$alpha
    }    
    sres
}

##' Improve MaxLik by Splitting
##'
##' This is a wrapper for \code{\link{estimateMaxLik}}, which see for
##' details about parameters and the return value
##' @title Splitting for Better BIC
##' @param V sticks
##' @param eta compositions
##' @param alpha tuning
##' @param params omega and psi
##' @param tab \code{wttab(...)}
##' @param ... other args to pass down
##' @param verbose print extra stuff, if 2 print lots
##' @return list with updated parms
##' @export
##' @author Charles Berry
splitMaxLik <- function(V,eta,alpha,params,tab,...,verbose=FALSE){
    result <- "split"
    save.res <- NULL
    while( result == "split" ){
	if (verbose) cat("s")
	res <- splitBestMaxLik(V,eta,alpha,params,tab,...)
	result <- res$action
	if (result=="split") {
	    eta <- res$res$eta
	    V <- res$res$V
	    alpha <- res$res$alpha
	    save.res <- res$res
	}
    }
    if (verbose) cat("\n")
    if (!is.null(save.res))
	save.res
    else
	estimateMaxLik(V,eta,alpha,params,tab,...)
}

##' Make Model Parsimonious by Merging for Best BIC
##'
##' This is a wrapper for \code{\link{estimateMaxLik}}, which see for
##' details about parameters and the return value
##' @title Merging to Better BIC
##' @param V sticks
##' @param eta compositions
##' @param alpha tuning
##' @param params omega and psi
##' @param tab \code{wttab(...)}
##' @param ... other args to pass down
##' @param verbose print extra stuff, if 2 print lots
##' @return list with updated parms
##' @export
##' @author Charles Berry
mergeMaxLik <- function(V,eta,alpha,params,tab,...,verbose=FALSE){
    result <- "merge"
    save.res <- NULL
    while( result == "merge" ){
	if (verbose) cat("m")
	res <- comboMaxLik(V,eta,alpha,params,tab,...)
	result <- res$action
	if (result=="merge") {
	    eta <- res$res$eta
	    V <- res$res$V
	    alpha <- res$res$alpha
	    save.res <- res$res
	}
    }
    if (verbose) cat("\n")
    if (!is.null(save.res))
	save.res
    else
	estimateMaxLik(V,eta,alpha,params,tab,...)
}


##' Do Best Merge (or not)
##'
##' This is a wrapper for \code{\link{estimateMaxLik}}, which see for
##' details about parameters and the return value
##' @title Merge One Pair for Better BIC
##' @param V sticks
##' @param eta compositions
##' @param alpha tuning
##' @param params omega and psi
##' @param tab \code{wttab(...)}
##' @param min.n smallest pseudo n for BIC
##' @param ... other args to pass down
##' @return list with updated parms
##' @importFrom utils combn
##' @export
##' @author Charles Berry
comboMaxLik <- function(V,eta,alpha,params,tab,min.n=2.0,...){
    T <- nrow(eta)
    eta.rows <- combn(T,2,simplify=FALSE)
    prob.z <- dZ.V(V)
    nparm <- ncol(eta)
    combo.res <-
	lapply(eta.rows,
	       function(x){
		   tab2 <- wttab.slice(x,V,eta,params,tab)
		   eta2 <- eta[x,]
		   V2 <- dV.Z(prop.table(prob.z[x]))
		   res2 <- estimateMaxLik(V2,eta2,alpha,params,tab2,...)
		   res1 <- estimateMaxLik(1.0,matrix(colMeans(eta2),nrow=1),
					  alpha,params,tab2,...)
		   bicdiff <- -2*( res2$logLik-res1$logLik ) +
		       log(max(min.n,sum(tab2$n)))*nparm
		   list(bic.diff=bicdiff,split=res2,merge=res1)
	       })
    worst.indx <- which.max(sapply(combo.res,"[[","bic.diff"))
    worst.bic <- combo.res[[worst.indx]]$bic.diff
    if (worst.bic<0){
	list(action="none",res=NULL)
    } else {
	merge.rows <- eta.rows[[worst.indx]]    
	eta.new <- rbind(eta[-merge.rows,],combo.res[[worst.indx]]$merge$eta)
	V.new <- dV.Z(c(prob.z[-merge.rows],sum(prob.z[merge.rows])))
	res.new <- estimateMaxLik(V.new,eta.new,alpha,params,tab,...)
	list(action="merge", res=res.new)
    }
}


##' Successively Delete Compositions Using AIC or Similar Criterion
##'
##' 
##' @title Delete Unneeded Compositions  
##' @param V sticks
##' @param eta compositions
##' @param alpha tuning
##' @param params omega and psi
##' @param tab \code{wttab(...)}
##' @param XIC multiplier. 2.0 yields AIC.
##' @param ... pass these also to \code{estimateMaxLik}
##' @param verbose if \code{TRUE} monitor progress
##' @return an object like \code{\link{estimateMaxLik}}
##' @importFrom parallel mclapply detectCores
##' @export
##' @author Charles Berry
deleteMaxLik <- function(V,eta,alpha,params,tab,XIC=2.0,...,verbose=FALSE){
    T <- nrow(eta)
    eta.rows <- 1:T
    nparm <- ncol(eta)
    res1 <- estimateMaxLik(V,eta, alpha,params,tab,...)
    action <- "merge"
    while (action=="merge"){
	combo.res <-
	    mclapply(eta.rows,
                     function(x){
                         eta2 <- res1$eta[-x,]
                         V2 <- dV.Z(prop.table(res1$prob.z[-x]))
                         res2 <- estimateMaxLik(V2,eta2,res1$alpha,params,tab,...)
                         bicdiff <- -2*( res1$logLik - res2$logLik ) +
                             XIC*nparm
                         list(bic.diff=bicdiff,merge=res2)
                     },mc.cores=round(detectCores()*3/4))
	worst.indx <- which.max(sapply(combo.res,"[[","bic.diff"))
	worst.bic <- combo.res[[worst.indx]]$bic.diff
	res.new <- combo.res[[worst.indx]][["merge"]]
	if (worst.bic<0){
	    action <- "none"
	    res.new <- res1
	} else {
	    if (verbose) cat("m")
	    res1 <- res.new
	    action <- if (length(res1$V)>2) "merge" else "none"
	    eta.rows <- 1:nrow(res1$eta)
	}
    }
    if (verbose) cat("\n")
    res.new
}





##' Do Best Split (or not)
##'
##' This is a wrapper for \code{\link{estimateMaxLik}}, which see for
##' details about parameters and the return value
##' @title Split Once for Better BIC
##' @param V sticks
##' @param eta compositions
##' @param alpha tuning
##' @param params omega and psi
##' @param tab \code{wttab(...)}
##' @param min.n smallest pseudo n for BIC
##' @param reflect Use random, antithetical compositions
##' @param ... other args to pass down
##' @return list with updated parms
##' @export
##' @author Charles Berry
splitBestMaxLik <- function(V,eta,alpha,params,tab,min.n=2.0,reflect=TRUE,...){
    T <- nrow(eta)
    eta.rows <- as.list(1:T)
    prob.z <- dZ.V(V)
    nparm <- ncol(eta)
    split.res <-
	lapply(eta.rows,
	       function(x){
		   tab2 <- wttab.slice(x,V,eta,params,tab)
		   eta1 <- eta[x,,drop=FALSE]
		   eta.trial <-
		       if (reflect)
		       {
			   eta.rnd <-
			       prop.table(pmax(rgamma(nparm,eta1*50),
					       .Machine$double.xmin)) 
			   eta.anti <- eta.reflect(eta.rnd,eta1)
			   rbind(eta.rnd,eta.anti)
		       } else {
			   reta(2)
		       }
		   res2 <- estimateMaxLik(dV.Z(len=2),eta.trial,
					  alpha,params,tab2,...)
		   res1 <- estimateMaxLik(1.0,eta1,alpha,params,tab2,...)
		   bicdiff <- -2*( res2$logLik-res1$logLik ) +
		       log(max(min.n,sum(tab2$n)))*nparm
		   list(bic.diff=bicdiff,split=res2,merge=res1)
	       })

    best.indx <- which.min(sapply(split.res,"[[","bic.diff"))
    best.bic <- split.res[[best.indx]]$bic.diff
    if (best.bic>0){
	list(action="none",res=NULL)
    } else {
	res.best <- split.res[[best.indx]]$split
	split.row <- eta.rows[[best.indx]]    
	eta.new <- rbind(eta[-split.row,],res.best$eta)
	V.new <- dV.Z(c(prob.z[-split.row],prob.z[split.row]*res.best$prob.z))
	res.new <- estimateMaxLik(V.new,eta.new,alpha,params,tab,...)
	list(action="split", res=res.new)
    }
}


expect.zw <- function(V, eta, params, tab){
    omega.psi <-
	params$omega %*% diag(params$psi)
    eodp <- prop.table(eta%*%omega.psi,1)
    prob.z <- dZ.V(V)
    like.zw <- t(dmulti(tab$tab,eodp))
    e.zw <- prop.table(like.zw*prob.z,2)%*%diag(tab$n)
}


wttab.slice <- function(rowindex,V,eta,params,tab,tol=1e-2){
    evals <- expect.zw(V, eta, params, tab)
    tab$n <-
	if (length(rowindex)==1)
	    evals[rowindex,]
	else
	    colSums(evals[rowindex,])
    tab$data.index <- NULL
    ok <- tab$n>tol
    if (!all(ok)) {
	tab$tab <- tab$tab[ok,]
	tab$n <- tab$n[ok]
    }
    tab}

eta.reflect <- function(x,base) prop.table(base/x)

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
    T <- nrow(eta)
    stopifnot(length(V)==T)
    omega <- params$omega
    lamb <- params$lambda
    s <- params$s
    stopifnot(length(s)==2)
    eps <- params$epsilon
    stopifnot(length(eps)==2)
    psi <- params$psi
    K <- ncol(tab$tab)
    stopifnot(length(psi)==K)
    stopifnot(length(lamb)==K)
    stopifnot(nrow(omega)==K)
    stopifnot(ncol(omega)==K)
    stopifnot(ncol(eta)==K)
    ka <- ko <- K
    ## `y' is used below where `W' might have better mnemonic
    ## meaning. This is a legacy of code written before subsampling was
    ## taken into account.
    W <- tab$tab
    Wplus <- rowSums(W)
    Nw <- length(Wplus)
    n <- tab$n
    ndat <- as.integer(nrow(W))
    stopifnot(ndat==length(n))
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
