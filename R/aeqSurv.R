#
# Create time values such that tiny differences are treated as a tie
#  The decision and tolerance are the same as all.equal
#
# see methods document: tied times
aeqSurv <- function(x, tolerance = sqrt(.Machine$double.eps)) {
    if (!missing(tolerance)) {
        if (!is.numeric(tolerance) || length(tolerance)!=1 ||
            !is.finite(tolerance))
            stop("invalid value for tolerance")
        if (tolerance <=0) return(x)  # do nothing
    }

    if (!is.Surv(x)) stop("argument is not a Surv object")
    y <- sort(unique(c(x[, -ncol(x)])))
    y <- y[is.finite(y)]  #someone may hand us an INF

    dy <- diff(y)
    tied <- ((dy <=tolerance) |( (dy/ mean(abs(y)) <=tolerance)))
    if (!any(tied)) return(x)   # all values are unique

    # There were ties.  Bin the data by the unique values that were found
    cuts <- y[c(TRUE, !tied)]  # set of unique values
    if (ncol(x) ==2) {  # simple survival
        z <- findInterval(x[,1], cuts)   # map each time point to an interval
        z <- cbind(cuts[z], as.integer(x[,2]))
    }
    else {
        z <- matrix(findInterval(x[,1:2], cuts), ncol=2)
        # We may have created zero length intervals
        zeros <- which(z[,1] == z[,2])
        if (length(zeros)>0 && any(x[zeros,1] != x[zeros,2]))
            stop("aeqSurv exception, an interval has effective length 0")
        z <- cbind(matrix(cuts[z], ncol=2), as.integer(x[,3]))
    }

    attributes(z) <- attributes(x)
    z
}
