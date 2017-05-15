##' Construct Lattice on a (hyper)Simplex 
##'
##' To mimic sampling from Dirichlet prior a lattice is made so that
##' samples from it approximate those from a Dirichlet( 1, 1, ... , 1)
##' distribution.
##' @title Lattice On A Simplex
##' @param from usually 0L - the origin of the lattice.
##' @param to integer giving the end of the last interval to be
##'     formed.
##' @param times the dimension of the simplex, i.e. \code{k} for a k-
##'     simplex.
##' @return matrix with \code{times+1} rows. The nubmer of columns is
##'     a function of the number of rows and of \code{to - from}.
##' @export
##' @author Charles Berry
##' @examples
##' simplexGrid(0, 3, 2)
simplexGrid <-
    function(from,to,times){
        times <- as.integer(times)
        stopifnot(times>=1,to>=from)
        if (times==1) {
            y <- 0:to
            rbind(y,to-y,deparse.level=0)
        } else {
            res <-
                sapply(from:to,
		       function(x)
                           rbind(x,
                                 simplexGrid(from=0,to=to-x,times=times-1L),
                                 deparse.level=0))
            if (is.list(res)) do.call(cbind,res) else res
        }
    }

##' Move Compositions to Simplex Lattice
##'
##' 
##' @title Nearest Lattice Point on a Simplex
##' @param p a matrix of compositions
##' @param incr numeric value st \code{0<incr && incr < 1}.  The
##'     inverse will be the number of intervals to use in forming a
##'     lattice on a simplex.
##' @return a matrix of compositions
##' @export
##' @author Charles Berry
roundSimplex <- function(p,incr=0.1)
{
    Nincr <- round(1/incr)
    start <- as.data.frame(p*Nincr)
    trial <- trunc(start)
    for (i in 1:ncol(p)){
        under <- round(Nincr-rowSums(trial))
        differ <- start-trial
        tops <- do.call(pmax,differ)
        trial <- trial +
            (under>0)*(differ==tops)
    }
    do.call(cbind,trial)/Nincr
}
