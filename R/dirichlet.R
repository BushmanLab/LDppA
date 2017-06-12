##' @importFrom stats rbinom rgamma

rdirichlet <-
    function(lam,log.q=FALSE)
{
    dm <- dim(lam)
    if (is.null(dm)) dm <- c(1,length(lam))
    res <- rgamma(length(lam),lam)
    dim(res) <- dm
    res <- quick.prop.table(res,1)
    if (log.q) log(res) else res
}

ddirichlet <-
    function(q,lam,log.p=FALSE)
{
    res <- 
	if (is.matrix(q)){
	    stopifnot(length(lam)==ncol(q))
	    (lam-1)%*%t(log(q)) +
		lgamma(sum(lam))-sum(lgamma(lam))
	} else {
	    stopifnot(length(lam)==length(q))
	    sum((lam-1)*log(q))+
		lgamma(sum(lam))-sum(lgamma(lam))
        }
    if (log.p) res else exp(res)
}

## use rmultinom instead??
rmulti <-
    function(n,p){
	stopifnot(length(n)==nrow(p))
	stopifnot(all(abs(rowSums(p)-1.0)<1e-12))
	res <- array(0L,dim(p))
	p0 <- 0
	for (i in 1:(ncol(p)-1)){
	    r <- rbinom(length(n),n,p[,i]/(1-p0))
	    res[,i] <- r
	    n <- n - r
	    p0 <- p0+p[,i]
	}
	res[,i+1] <- n
	res
    }

dmulti <-
    function(y,p,min.p=0.0){
	## for each row of y and row of p compute prob
	## result is nrow(y) X ncol(p) 
	if (!is.matrix(y)) dim(y) <- c(1,length(y))
	if (!is.matrix(p)) dim(p) <- c(1,length(p))
	nry <- nrow(y)
	nrp <- nrow(p)
	stopifnot(abs(rowSums(p)-1.0)<sqrt(.Machine$double.eps))
	res <- tcrossprod(y, log(p))
	res.M <- lgamma(rowSums(y)+1)-rowSums(lgamma(y+1))
	if (min.p==0.0) 
	    exp(res.M+res)
	else
	    pmax( exp(res.M+res) , min.p)
    }
