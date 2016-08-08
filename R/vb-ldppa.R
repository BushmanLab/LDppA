##' Variational Bayes for Latent Dirichlet Process Prior Allocation
##'
##' LDDPA taking alpha as a latent variable
##' @title ldppa.vb
##' @param wttab a \code{list} with elements \code{tab} a table of
##'     counts and \code{n} rowwise case weights for the table.
##' @param params original parameter values
##' @param initialValues intial values for variation parameters
##' @return \code{list} of updated variational parameters
##' @export
##' @author Charles Berry
ldppa.vb <- function(wttab,params,initialValues){
    tab <- wttab$tab
    n <- wttab$n
    s <- params$s
    omega <- params$omega
    lambda <- params$lambda
    K <- ncol(omega)
    N <- nrow(tab)
    w <- initialValues$w
    alpha <- w[1]/w[2]
    gamma <- initialValues$gamma
    tau <- initialValues$tau
    phi <- initialValues$phi
    Ntau <- nrow(initialValues$tau) # much smaller than N
    EqlogEta <- phi %*% ( digamma(tau) - digamma(rowSums(tau)) ) # [n,r]
    mutab <- rep(omega,each=N) * exp(EqlogEta)[,rep(1:K,K)] # [n,r,s]
    dim(mutab) <- c(N,K,K)
    mutab <- quick.prop.table(mutab,c(1,3))
    gamma <- cbind(1+colSums(phi*n)[-Ntau],
                   alpha + sum(phi*n) - cumsum(colSums(phi*n)[-Ntau]))
    w <- c(s[1]+Ntau-1, s[2]-sum(digamma(gamma[,2])-digamma(rowSums(gamma))))
    tau <- rep(lambda,each=Ntau) + t(phi)%*%(x.n.plus.r.plus(tab,mutab)*n)
    phi <- phi.initial(tab,mutab,gamma,tau)
    list(gamma=gamma,tau=tau,phi=phi,mu=mutab,w=w,alpha=alpha)
}

##' Fit Variational Bayes for Latent Dirichlet Process Prior Allocation
##'
##' This is a wrapper for \code{ldppa.vb}
##' @title lddpa.VarBayes
##' @param tab same as for \code{ldppa.vb}
##' @param params same as for \code{ldppa.vb}
##' @param initialVals an integer giving the number of mixture
##'     components or a list of initial values. See \code{ldppa.vb}.
##' @param n.iter an integer for the number of iterations (less 1)
##' @return the value of \code{ldppa.vb} plus an component giving the
##'     fit of the first interation
##' @export
##' @author Charles Berry
ldppa.VarBayes <-
    function(tab,params,initialVals=100L,n.iter=100L)
{
    ## omega,lambda,s in params
    if (is.numeric(initialVals)){
        omega <- params[["omega"]]
        lambda <- params[["lambda"]]
        s <- params[["s"]]
        TEE <- initialVals
        gam0 <- gamma.initial(TEE)
        tau0 <- tau.initial(tab,omega,lambda,TEE)
        initialVals <- list(
            gamma = gam0,
            w = cbind(s[1]+TEE-1,
  		    s[2] - sum(digamma(gam0[,2]) -
                                 digamma(rowSums(gam0)))),
            tau = tau0,
            phi = phi.initial(tab$tab,
  			    mu.initial(tab,omega) , gam0, tau0 ))
    }

    fit2 <- fit1 <- ldppa.vb(tab, params,initialVals)
    for (i.iter in 1:n.iter) fit2 <- ldppa.vb(tab, params, fit2)
    c(fit2,list(start=fit1))
}

##' Latent Dirichlet Process Prior Variational Bound on Log Posterior
##'
##' This is a wrapper for \code{lb}
##' @title lower.bound
##' @param tab as per \code{ldppa.vb}
##' @param params as per \code{ldppa.vb}
##' @param fit the result of \code{ldppa.vb}
##' @return The log posterior as determined by the variational
##'     approximation.
##' @export
##' @author Charles Berry
lower.bound <- function(tab,params,fit)

{
    lb(tab,
        fit[["w"]],
        fit[["gamma"]],
        fit[["tau"]],
        fit[["phi"]],
        fit[["mu"]],
        params[["lambda"]],
        params[["omega"]],
        fit[["alpha"]],
        params[["s"]])
}
