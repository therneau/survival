#
# After 20 years of use, some deficiencies in the layout of survfit objects
#  have become obvious, in particular that the starting point of time=0 
#  survival =1 is not included in the curves.  This has led to extra arguments
#  like start.time and p0. Hindsight is 20/20, as they say, and I can now see
#  that leaving that out was a mistake.
# This routine converts to a newer design called version 3 which inserts those
#  time points. It makes downstream processing of the object easier, and is used
#  by the plot functions, for instance.
# Whether survival version 3 will use this form natively is still not decided.
#
survfit23 <- function(x) {
    if (!inherits(x, "survfit")) stop("function requires a survfit object")
    if (!is.null(x$version3) && x$version3) return(x)  # already in 3.x format
    if (is.null(x$start.time)) start.time <- 0 else start.time <- x$start.time

    if (is.null(x$strata)) insert <- 1   # where to add the zero
    else insert <- unname(1 + cumsum(c(0, x$strata[-length(x$strata)])))
    
    same <- x$time[insert] == start.time  # no row need to be inserted here
    insert <- insert[!same]
    if (length(insert)==0) {  # nothing much to do
        drop <- c("start.time", "p0")
        newx <- unclass(x)[is.na(match(names(x), drop))]
        newx$version3 <- TRUE
        class(newx) <- class(x)
        return(newx)
    }
    if (!is.null(x$strata)) {
        newstrat <- x$strata
        newstrat[!same] <- newstrat[!same] +1
    }

    # actual work needs to be done
    # this adds rows to a vector or matrix, and preserves double/integer
    addto <- function(x, i, z) {
        # i = where to add, z = what to add
        n.add <- length(i)
        i2 <- i + 1:n.add -1

        if (is.matrix(x)) {
            # indx is the new rows that are equal to the old ones
            indx <- seq(1, n.add + nrow(x))[ -i2]
            newx <- matrix(x[1], nrow(x) + n.add, ncol(x))
            newx[indx,] <- x
            newx[i2,] <- z
        } else { 
           indx <- seq(1, n.add + length(x))[ -i2]
           newx <- rep(x[1], length(x) + n.add)
           newx[indx] <- x
           newx[i2] <- z
        }   
        newx
    }
    
    # I want the result to be a list in the same order
    newname <- names(x)[!(names(x) %in% c("start.time", "p0"))]
    new <- vector("list", length(newname))
    names(new) <- newname

    add1 <- c("surv", "lower","upper")
    add0 <- c("n.event", "n.censor", "n.add", "std.err", "cumhaz", "std.chaz")
    for (i in names(new)) {
        if (i=="time") new[[i]] <- addto(x[[i]], insert, start.time)
        else if (i=="n.risk") new[[i]] <- addto(x[[i]], insert, x$n.risk[1])
        else if (i=="pstate") new[[i]] <- addto(x[[i]], insert, x$p0)
        else if (i=="strata") new[[i]] <- newstrat
        else if (i %in% add0) new[[i]] <- addto(x[[i]], insert, 0L)
        else if (i %in% add1) new[[i]] <- addto(x[[i]], insert, 1L)
        else new[[i]] <- x[[i]]
    }

    new$version3 <- TRUE
    class(new) <- class(x)
    new
}
            
        
    
