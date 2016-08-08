##' @importFrom stats runif
scale.dirichlet <- function(x,scale) (x-1)/scale+1  

##' Latent Variables are sampled from their variational distributions
##' according to fitted values of variational parameters.
##'
##' 
##' @title Sample from Variational Distribution Fit
##' @param fit the result of \code{\link{ldppa.VarBayes}}
##' @param scale numeric value.  Values greater than 1.0 will tend to
##'     result in higher dispersion of the sampled values.  Ignored
##'     for \code{rQalpha}.
##' @return a single sample from the distribution
##' @export
##' @author Charles Berry
rQeta <- function(fit,scale=1) rdirichlet(scale.dirichlet(fit$tau,scale))

rQz <- function(fit,tab)
{
    di <- tab$data.index
    rowSums(t(apply(fit$phi,1,cumsum)[,di])<runif(length(di)))+1
}

##' @rdname rQeta
##' @export
rQV <- function(fit,scale=1)
{
    gam <- scale.dirichlet(fit$gamma,scale)
    c(rbeta(nrow(gam),gam[,1],gam[,2]),1)
}

##' @rdname rQeta
##' @export
rQalpha <- function(fit, scale=NULL)
{
    rgamma(1,fit$w[1],rate=fit$w[2])
}
rQx <- function(fit,tab)
{
    tabl <- tab$tab
    res <-
        sapply(tab$data.index,
     	  function(x)
     	  {
                   sapply(1:ncol(tabl),
                          function(y) if (tabl[x,y]==0)
                                          rep(0,ncol(tabl))
     				 else
                                          rmulti(tabl[x,y],fit$mu[x,,y]))
     	  })
    dim(res) <- c(ncol(tabl),ncol(tabl),length(tab$data.index))
    aperm(res,c(3,1,2))
}
