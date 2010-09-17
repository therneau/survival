# $Id: survreg.control.S 11236 2009-02-14 11:46:53Z therneau $
survreg.control <- function(maxiter=30, rel.tolerance=1e-9, 
			    toler.chol=1e-10, iter.max, debug=0,
			    outer.max = 10) {

    if (missing(iter.max)) {
	iter.max <- maxiter
	}
    else  maxiter <- iter.max
    list(iter.max = iter.max, rel.tolerance = rel.tolerance, 
	 toler.chol= toler.chol, debug=debug,
	 maxiter=maxiter, outer.max=outer.max)
    }
