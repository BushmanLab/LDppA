##' Initialize Gibbs Sampler
##'
##' Use approximate maximum a posteriori estimates to initialize a
##' Gibbs sampler.  Find the maximum likelihood values, cluster a
##' sample of the values to find eta and V, then get alpha.
##' @title eva.initialize
##' @param tab a \code{list} with elements \code{tab} a table of
##'     counts and \code{n} rowwise case weights for the table.
##' @param omega the confusion matrix possibly scaled by the
##'     subsampling fractions.
##' @param lambda Dirichlet parameter, usually the unit vector .
##' @param k number of rows eta should contain
##' @param max.rows hard limit on number of rows to sample.
##' @param cap.wt whether to truncate the number of single cell ISs to the
##'     largest number of 2 cell ISs in sampling combinations
##' @return list with elements \code{eta}, \code{V}, and \code{alpha}.
##' @export
##' @author Charles Berry
eva.initialize <-
    function(tab,omega,lambda,k=NULL,max.rows=NULL, cap.wt=TRUE) {
        eta.from.phi <-
            function(phi) c(1,exp(phi))/(1+sum(exp(phi)))
        argmax.llk <- function(phi,w,omega.psi,tol=1e-10){
            log.eta <- c(0,phi)-max(0,phi)
            eta <- exp(log.eta)/sum(exp(log.eta))
            eta.op <- eta%*%omega.psi
            eta.op <- eta.op/sum(eta.op)
            res <- dmultinom(w,prob=eta.op,log=TRUE)
            res
        }
        y <- tab$tab
        wt <- tab$n
        cells <- rowSums(y)
        boot.wt <-
            if (cap.wt) pmin(wt, max(wt[cells>1])) else wt
        if (is.null(max.rows)) max.rows <-
                                   min(sum(wt),k*100)
        if (is.null(k)) k <- length(lambda)
        boot.samp <- sample(nrow(y),max.rows,replace=TRUE,prob=boot.wt)
        rowvals <- apply(y,1,
                         function(x){
                             opt <- optim(rep(0,length(lambda)-1),
                                          argmax.llk,
                                          w=x,
                                          omega=omega,
                                          control=list(fnscale=-1))
                             eta.from.phi(opt$par)
                         })
        y <- t(rowvals[,boot.samp])*(wt*cells)[boot.samp]
        ind <- cutree(hclust(dist(y)),k=k)
        tau <- unname(coef(lm(y~0+as.factor(ind)))+rep(lambda,each=k))
        eta <- prop.table(tau,1)
        pr.z <- prop.table(rowSums(tau))
        v <- pr.z/rev(cumsum(rev(pr.z)))
        v[length(v)] <- 1.0
        alpha <- -(k-1)/sum(log(1-head(v,-1)))
        list(alpha=alpha,eta=eta,V=v)
    }
