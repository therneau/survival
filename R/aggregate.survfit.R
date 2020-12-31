#
# Aggregate function for survfit objects, used to get population
#  averages
#
aggregate.survfit <- function(x, by=NULL, FUN= mean, ...) {
    if (!inherits(x, "survfit")) stop("x must be a survfit object")
    
    dd <- dim(x)
    dd.data <- dd["data"]
    if (is.null(dd.data)) stop("survfit object does not have a 'data' margin")
    if (is.null(by)) index <- rep.int(1L, dd.data)
    else {
        if (is.list(by)) {
            blen <- sapply(by, length)
            if (any(blen != dd.data)) stop("arguments must have the same length")
        }
        else {
            if (length(by)!= dd.data) stop("arguments must have the same length")
            by <- list(by)
        }

        # create an integer index that will apply to each column of surv, pstate,
        #  or cumhaz
        index <- tapply(by[[1]], by)    # integer version of "by"
        index <- match(index, sort(unique(index))) # no holes in the sequence
        if (all(index == index[1])) by <- NULL  # all in one group
    }

    # test that FUN is okay, using a dummy vector of the right length
    test <- tapply(seq.int(dd.data), index, FUN)
    if (is.list(test) || length(test) != max(index) || !is.numeric(test))
        stop("FUN must return a single value summary")

    # these components don't collapse
    j <- match(c("std.err", "std.cumhaz", "lower", "upper", "conf.int",
                       "conf.type", "logse", "cumhaz"), names(x), nomatch= 0)
    newx <- unclass(x)[-j]

    if (is.null(by)) { # simple case
        if (!is.null(x$surv)) {
            if (missing(FUN)) newx$surv <- rowMeans(x$surv)
            else              news$surv <- apply(x$surv, 1, FUN)
        }
        if (!is.null(x$pstate)) 
            newx$pstate <- apply(x$pstate, c(1,3), FUN)
    }
    else {
        if (FALSE) {
        #if missing(FUN)) { # use a fast algorithm tailored to the mean
            if (!is.null(x$surv))
                newx$surv <- .Call(Csurvfitmean, x$surv, dim(x$surv), index- 1L)
            if (!is.null(x$pstate))
                newx$pstate <-.Call(Csurvfitmean, x$pstate, dim(x$pstate), 
                                    index -1L)
        }
        else {  # the complicated one
            if (!is.null(x$surv)) {
                temp <- apply(x$surv, 1, function(z) tapply(z, by, FUN))
                newx$surv <- t(temp)
            }
            if (!is.null(x$pstate)) {
                temp <- apply(x$pstate, c(1,3), function(z) tapply(z, by, FUN))
                newx$pstate <- aperm(temp, c(2,1,3))
            }
        }
    }    

    if (is.null(by)) newx$newdata <- NULL
    else newx$newdata <- data.frame("(aggregate)" = seq.int(max(index)))

    class(newx) <- class(x)
    newx
}
