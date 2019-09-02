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
    if (!is.null(x$version) && x$version==3) return(x)  # already in 3.x format
    if (is.null(x$start.time)) start.time <- 0 else start.time <- x$start.time

    if (is.null(x$strata)) insert <- 1   # where to add the zero
    else insert <- unname(1 + cumsum(c(0, x$strata[-length(x$strata)])))
    
    same <- x$time[insert] == start.time  # no row need to be inserted here
    insert <- insert[!same]
    if (length(insert)==0) {  # nothing much to do
        drop <- c("start.time", "p0")
        new <- unclass(x)[is.na(match(names(x), drop))]
        new$version <- 3
        class(new) <- class(x)
        return(new)
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
        }
        else{ 
           indx <- seq(1, n.add + length(x))[ -i2]
           newx <- rep(x[1], length(x) + n.add)
           newx[indx] <- x
           newx[i2] <- z
        }   
        newx
    }
    
    # I want the result to be a list in the same order
    newname <- names(x)[!(names(x) %in% c("start.time", "p0", "sp0"))]
    new <- vector("list", length(newname))
    names(new) <- newname

    add1 <- c("surv", "lower","upper")
    add0 <- c("n.event", "n.censor", "n.add", "cumhaz", "std.chaz")

   if (is.null(x$sp0)) sp0 <- 0 else sp0 <- x$sp0
    if (!is.null(x$p0)) {
      if (any(same)) {# we have to subscript p0 and sp0
            # if p0 isn't a matrix, we can't end up here BTW
            p00 <- x$p0[!same,]
            if (!is.null(x$sp0)) sp0 <- x$sp0[!same,]
      }
      else p00 <- x$p0
    }

    for (i in names(new)) {
        if (i=="time") new[[i]] <- addto(x[[i]], insert, start.time)
        else if (i=="n.risk") {
            if (is.matrix(x$n.risk))
                new[[i]] <- addto(x[[i]], insert, x$n.risk[insert,])
            else new[[i]] <- addto(x[[i]], insert, x$n.risk[insert])
        }
        else if (i=="pstate") new[[i]] <- addto(x[[i]], insert, p00)
        else if (i=="strata") new[[i]] <- newstrat
        else if (i=="std.err") new[[i]] <- addto(x[[i]], insert, sp0)
        else if (i %in% add0) new[[i]] <- addto(x[[i]], insert, 0L)
        else if (i %in% add1) new[[i]] <- addto(x[[i]], insert, 1L)
        else new[[i]] <- x[[i]]
    }
    
    if (!inherits(x, "survfitms")) {
        addcol <- function(x) cbind(0, x)
        # we need to fix up influence objects, which have subjects as the
        #  rows and time as the columns.  If there are multiple curves it
        #  will be a list with one element per curve
        for (i in c("influence.surv", "influence.chaz")) {
            if (!is.null(x[[i]])) {
                if (is.list(x[[i]])) new[[i]] <- lapply(x[[i]], addcol)
                else new[[i]] <- addcol(x[[i]])
            }
        }
    }
                     
    if (is.null(new$logse)) {
        # reprise the logic of the older code
        if (inherits(x, "survfitms")) x$logse <- FALSE
        else x$logse <- TRUE
    }
    
    if (is.null(x$cumhaz) && class(x)[1] == "survfit") {  # fill it in!
        new$cumhaz <- -log(new$surv)
        if (!is.null(x$std.err)) {
            if (x$logse) new$std.chaz <- new$std.err
            else         new$std.chaz <- new$std.err/new$surv
        }
    }
    new$version <- 3
    class(new) <- class(x)
    new
}
            
survfit32 <- function(x) {
    if (!inherits(x, "survfit")) stop("function requires a survfit object")
    if (is.null(x$version) || x$version<3) return(x)  # already in proper form
    
    if (is.null(x$strata)) first <- 1
    else {
        last <- cumsum(x$strata)
        first <- 1+ c(0, last[-length(last)])
        x$strata <- x$strata -1L
    }

    x$start.time <- x$time[1]
    for (i in c("time", "n.risk", "n.event", "n.censor", "cumhaz",
                "std.chaz", "lower", "upper")){
        if (!is.null(x[[i]])) {
            if (is.matrix(x[[i]])) x[[i]] <- x[[i]][-first,,drop=FALSE]
            else x[[i]] <- x[[i]][-first]
        }
    }       

    if (inherits(x, "survfitms")) {
        if (is.matrix(x$pstate)) {
            x$p0 <- x$pstate[first,]
            x$pstate <- x$pstate[-first,,drop=FALSE]
            if (!is.null(x$std.err)) {
                x$sp0 <- x$std.err[first,]
                x$std.err <- x$std.err[-first,, drop=FALSE]
            }
        }else {
            x$p0 <- x$pstate[first]
            x$pstate <- x$pstate[-first]
            if (!is.null(x$std.err)) {
                x$sp0 <- x$std.err[first]
                x$std.err <- x$std.err[-first]
            }
        }   
    } else {
        if (is.matrix(x$surv)) {
            x$surv <- x$surv[-first,,drop=FALSE]
            if (!is.null(x$std.err)) {
                x$std.err <- x$std.err[-first,, drop=FALSE]
            }
        }else {
            x$surv <- x$surv[-first]
            if (!is.null(x$std.err)) {
                x$std.err <- x$std.err[-first]
            }
        }   
    } 

    x$version <- 2
    x
}

    
