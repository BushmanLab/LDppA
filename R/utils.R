prob.z.v <- function(v) v*exp(c(0,cumsum(log1p(-v)[-length(v)])))

logP.ztab.v <- function(ztab,v){
    pie <- log(prob.z.v(v))
    sum(pie*ztab)
}
