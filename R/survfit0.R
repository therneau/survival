# This is a utility function, used by printing and plotting.
#   If the original call had time0=FALSE as the option, then the curve does not
# have time 0 in the surv, pstate or cumhaz, and corresponding std dev matrices
# This function's job is to turn it into the time0=TRUE form.
# For ordinary survival surv=1, cumulative hazard =0, and standard errors
# of 0.  For a multi-state curve the values of p0 and sd0 are used to fill
# in these values.  
#
#  The influence matrix, if present, is also filled out using i0, if present.
# Older versions of survfit had included time0 in the influence.
#
survfit0 <- function(x, ...) {
    if (!inherits(x, "survfit")) stop("function requires a survfit object")
    # the !is.na is for backwards compatability with a saved object
    if (inherits (x, "survfit0") || (!is.null(x$time0) && x$time0)) return(x)

    t0 <- x$t0 # present for multi-state
    if (is.null(t0)) t0 <- min(c(0, x$time))  # ordinary survival
 
    # insert will point to the first row of each stratum
    if (is.null(x$strata)) insert <- 1   # where to add the zero
    else insert <- unname(1 + cumsum(c(0, x$strata[-length(x$strata)])))

    same <- (x$time[insert] == t0)  # no row need to be inserted here
    insert <- insert[!same]
    if (length(insert)==0) return(x)   # nothing  to do
    # actual work needs to be done
 
    if (!is.null(x$strata)) {
        newstrat <- x$strata
        newstrat[!same] <- newstrat[!same] +1
    }

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
            newx[i2,] <- rep(z, each=dd[2])
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
    # change on 12/2020: p0 and sp0 will be redundant, but leave them in
    #  anyway.  It makes it easier for survfitCI and residuals.survfit.
    # Note on 3/2024: actually, p0 is not redundant when there is an event
    #  exactly at the starting time, but sp0 is not needed
    newname <- names(x)[!(names(x) %in% c("sp0"))]
    new <- vector("list", length(newname))
    names(new) <- newname

    add1 <- "surv"
    add0 <- c("n.event", "n.censor", "n.enter", "n.transition",
              "cumhaz", "std.chaz", "std.auc")
    if (inherits(x, "survfitms")) add0 <- c(add0, "lower", "upper")
    else add1 <- c(add1, "lower", "upper")

    if (!is.null(x$p0)) {  # multi-state
        if (is.null(x$sp0)) sp0 <- 0 else sp0 <- x$sp0
        if (any(same)) {# we have to subscript p0 and sp0
            # if p0 isn't a matrix, we can't end up here BTW
            p00 <- x$p0[!same,,drop=FALSE]
            if (!is.null(x$sp0)) sp0 <- x$sp0[!same,,drop=FALSE]
        }
        else p00 <- x$p0
    }
    else { # this call is simply adding a time 0 to the front
        if (!is.null(x$xtd.err)) sp0 <- x$std.err[insert,]
        else sp0 <- 0
    }

    for (i in names(new)) {
        if (i=="time") new[[i]] <- addto(x[[i]], insert, t0)
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
    } else if (!is.null(x$influence)) { # this is harder, multi-state
        addi0 <- function(x, i0 =0) {
            # the influence matrix will be (observations, times, states),
            # we need to add a new time
            dd <- dim(x)
            temp <- array(0., dim=c(dd[1], 1L + dd[2], dd[3]))
            temp[,-1,] <- x
            if (!is.null(i0)) temp[,1,] <- i0
            temp
        }
        if (is.list(x$influence)) {
            if (is.null(x$i0)) new$influence <- lapply(x$influence, addi0)
            else { # i0 will be a list as well
                new$influence <- vector("list", length(x$influence))
                for (i in 1:length(x$influence))
                    new$influence[[i]] <- addi0(x$influence[[i]], x$i0[[i]])
            }
        } else new$influence <- addi0(x$influence, x$i0)
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
