# This is a utility function, used by printing and plotting.
#    It's job is to add a "time 0" to the survival curve, and at the same time
# to fill in values for the responses at that time.  For ordinary survival
# the standard response is surv=1, cumulative hazard =0, and standard errors
# of 0.  For a multi-state curve the values of p0 and sd0 are used to fill
# in these values.  
#   The influence matrix, if present, is also filled out to have a 0's for the
# time, or for a multi-state object, use the influence matrix as provided.
#
# It is legal for start.time to be a vector, so that multiple curves start in
#  different places.  I don't yet have a use for that, but had considered one.
#
survfit0 <- function(x, start.time =0) {
    if (!inherits(x, "survfit")) stop("function requires a survfit object")
    if (inherits (x, "survfit0")) return(x)

    start.time <- min(start.time, 0, x$time)  #start.time might be NULL

    if (is.null(x$strata)) insert <- 1   # where to add the zero
    else insert <- unname(1 + cumsum(c(0, x$strata[-length(x$strata)])))

    same <- x$time[insert] == start.time  # no row need to be inserted here
    insert <- insert[!same]
    if (length(insert)==0) return(x)   # nothing  to do

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
            if (length(colnames(x)) >0) colnames(newx) <- colnames(x)
            newx[indx,] <- x
            newx[i2,] <- z
        } else if (is.array(x)) {
            # pstate is sometimes an array
            dd <- dim(x)
            indx <- seq(1, n.add + dd[1])[ -i2]
            newx <- matrix(x[1], dd[1] + n.add, dd[2]*dd[3])
            newx[indx,] <- c(x)
            newx[i2,] <- z
            dim(newx) <- c(dd[1]+n.add, dd[2], dd[3])
            dimnames(newx) <- list(NULL, NULL, dimnames(x)[[3]])
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

    if (!is.null(x$p0)) {
        if (is.null(x$sp0)) sp0 <- 0 else sp0 <- x$sp0
        if (any(same)) {# we have to subscript p0 and sp0
            # if p0 isn't a matrix, we can't end up here BTW
            p00 <- x$p0[!same,]
            if (!is.null(x$sp0)) sp0 <- x$sp0[!same,]
        }
        else p00 <- x$p0
    }
    else { # this call is simply adding a time 0 to the front
        p00 <- x$pstate[insert,]
        if (!is.null(x$xtd.err)) sp0 <- x$std.err[insert,]
        else sp0 <- 0
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
    class(new) <- c("survfit0", class(x))
    new
}
