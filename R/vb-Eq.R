## omit the unused Tth values (for which EqlogV=0 and Eqlog1mV=-Inf)
EqlogV <- function(gamma) digamma(gamma[,1])-digamma(rowSums(gamma))

Eqlog1mV <- function(gamma) digamma(gamma[,2])-digamma(rowSums(gamma))


EqlogPAlphamqAlpha <- function(s,w){
    (s[1]-w[1])*(digamma(w[1])-log(w[2]))-
        (s[2]-w[2])*w[1]/w[2] +
        s[1]*log(s[2])-w[1]*log(w[2]) -
        lgamma(s[1]) + lgamma(w[1])}


EqlogPVmqV <- function(gamma,alpha)
    sum(
        log(alpha)+(alpha-1)*Eqlog1mV(gamma),
        -(lgamma(rowSums(gamma))-rowSums(lgamma(gamma))+
          (gamma[,1]-1)*EqlogV(gamma) +
          (gamma[,2]-1)*Eqlog1mV(gamma)))


## rows of tau contain dirichlet params
EqlogPEtamqEta <- function(tau,lambda){
    sum(lgamma(sum(lambda))-sum(lgamma(lambda))+
        (lambda-1)%*%t(digamma(tau)-digamma(rowSums(tau)))+
        ( -lgamma(rowSums(tau))+rowSums(lgamma(tau))-
          rowSums((tau-1)*(digamma(tau)-digamma(rowSums(tau))))
        ))
}

EqlogPZmqZ <- function(tab,phi,gamma){
    n <- tab$n
    t.phiplus <- apply(phi,1,function(x) sum(x)-cumsum(x))
    sum(c(Eqlog1mV(gamma),0)%*%t.phiplus%*%n)+
        sum(n%*%phi%*%c(EqlogV(gamma),0))-
        sum((phi*log(phi)*n)[phi>0])
}

EqlogPXmqX <- function(tab,mu,phi,tau){
    y <- tab$tab
    n <- tab$n
    K <- ncol(tau)
    N <- nrow(phi)
    ElogEta <- phi%*%(digamma(tau)-digamma(rowSums(tau))) # N x K
    mu.y <- mu*as.vector(y[,rep(1:K,each=K)]*n) # N x K x K
    sum(rowSums(mu.y,dims=2)*ElogEta) -sum(mu.y*log(mu))
}



###


EqlogPY <- function(tab,mu,omega){
    y <- tab$tab
    n <- tab$n
    N <- nrow(y)
    K <- ncol(y)
    sum(as.vector(y[,rep(1:K,each=K)]*n)*mu*rep(log(omega),each=N))
}


lb <- function(tab,w,gamma,tau,phi,mu,lambda,omega,alpha,s){
    EqlogPAlphamqAlpha(s,w)+ 
        EqlogPY(tab,mu,omega)+
        EqlogPXmqX(tab,mu,phi,tau)+
        EqlogPZmqZ(tab,phi,gamma)+
        EqlogPEtamqEta(tau,lambda)+
        EqlogPVmqV(gamma,alpha)
}
