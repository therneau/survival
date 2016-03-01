# Create the Aalen-Johansen estimate by joining a set of 
#  survival curves.  (Actually the cumulative hazard estimates are used.)
#
survfit.matrix <- function(curves, p0) {
    Call <- match.call()
    if (!is.matrix(curves)) 
             stop("input must be a square matrix of survival curves")
    if (!is.list(curves)) stop("input must be a matrix of survival curves")
    if (nrow(curves) != ncol(curves)) 
            stop("input must be a square matrix survival curves")
    nstate <- nrow(curves)

    nonzero <- !sapply(curves, is.null) # no transition here
    if (any(sapply(curves[nonzero], function(z) !inherits(z, 'survfit'))))
            stop("input must be a square matrix survival curves")
    if (sum(nonzero) < 2)
        stop("input must have at least 2 transitions")
    
    classes <- lapply(curves[nonzero], class)
    # Make sure we were sent the right things.  If any of the curves inherit
    #  from "survfitms", they are the result of a prior AJ computation;
    #  such recursion does not lead to valid estimates.
    if (any(sapply(classes, function(x) inherits(x, "survfitms"))))
        stop("multi-state curves are not a valid input")
    if (all(sapply(classes, function(x) length(x)==1 && x=='survfit')))
        type <- 'survfit'
    else if (all(sapply(classes, function(x) inherits(x, 'survfit.coxph'))))
        type <- 'survfit.coxph'
    else stop ("all curves must be of the same type")

    utime <- lapply(curves[nonzero], function(x) x$time[x$n.event>0])
    utime <- sort(unique(utime))  # set of unique times
    cumhaz<- lapply(curves[nonzero], function(x)
                    summary(x, times=utime, extend=TRUE)$cumhaz)
    jumps <- matrix(unlist(lapply(cumhaz, function(x) diff(c(0, x)))),
                    ncol= sum(nonzero))

    Tmat <- diag(nstate)
    if (missing(p0)) p0 <- c(1, rep(0, nstate-1))
    else if (sum(p0) !=0 || any(p0<1) || length(p0) != nstate)
        stop("invalid p0 vector")

    prev  <- matrix(0., nrow= 1+length(utime), ncol=nstate)
    prev[1,] <- p0
    for (i in 1:length(utime)) {
        Tmat[nonzero] <- jumps[i,]
         prev[i+1,] <- prev[i,] %*% Tmat
    }
    
    list(time=c(0, utime), P = prev)
}
#
