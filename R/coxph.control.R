#
# Gather all of the control parameters for coxph into one spot
#
coxph.control <- function(eps=1e-9, 
                          toler.chol = .Machine$double.eps ^ .75, 
			  iter.max=20,
			  toler.inf= sqrt(eps), outer.max=10,
                          timefix =TRUE) {
    if (!is.numeric(iter.max) ||iter.max <0) stop("Invalid value for iterations")
    if (!is.numeric(eps) || eps <=0) stop ("Invalid convergence criteria")
    if (!is.numeric(toler.chol) || toler.chol <=0)
        stop("invalid value for toler.chol")
    if (!is.numeric(eps) || eps <=0) stop("eps must be > 0")
    if (eps <= toler.chol) 
	    warning("For numerical accuracy, tolerance should be < eps")
    if (!is.numeric(toler.inf) || toler.inf <=0) 
        stop ("The toler.inf setting must be >0")
    if (!is.logical(timefix)) stop("timefix must be TRUE or FALSE")
    if (!is.numeric(outer.max) || outer.max <=0) 
        stop("invalid value for outer.max")
    list(eps=eps, toler.chol=toler.chol, iter.max=as.integer(iter.max), 
	 toler.inf=toler.inf, outer.max=as.integer(outer.max), 
         timefix=timefix)
    }
