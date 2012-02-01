#
# Gather all of the control parameters for coxph into one spot
#
coxph.control <- function(eps=1e-9, 
                          toler.chol = .Machine$double.eps ^ .75, 
			  iter.max=20,
			  toler.inf= sqrt(eps), outer.max=10 ) {
    if (iter.max <0) stop("Invalid value for iterations")
    if (eps <=0) stop ("Invalid convergence criteria")
    if (eps <= toler.chol) 
	    warning("For numerical accuracy, tolerance should be < eps")
    if (toler.inf <=0) stop ("The inf.warn setting must be >0")
    list(eps=eps, toler.chol=toler.chol, iter.max=as.integer(iter.max), 
	 toler.inf=toler.inf, outer.max=as.integer(outer.max))
    }
