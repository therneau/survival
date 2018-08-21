# Create the Aalen-Johansen estimate by joining a set of 
#  survival curves.  Actually the cumulative hazard estimates are used.
#
survfit.matrix <- function(formula, p0, method=c("discrete", "matexp"), 
                           start.time, ...) {
    Call <- match.call()
    curves <- formula
    if (!is.matrix(curves)) 
             stop("input must be a square matrix of survival curves")
    if (!is.list(curves)) stop("input must be a matrix of survival curves")
    if (nrow(curves) != ncol(curves)) 
            stop("input must be a square matrix survival curves")
    nstate <- nrow(curves)

    states <- names(p0)
    if (is.null(states)) {
        dname <- dimnames(curves)
        if (!is.null(dname)[[1]]) states <- dname[[1]]
        else if (!is.null(dname[[2]])) states <- dname[[2]]
        else states <- 1:nstate
    }

    nonzero <- sapply(curves, function(x) length(x) > 0) # transitions
    curves <- curves[nonzero]  #toss away the NULLS
    if (any(sapply(curves, function(z) !inherits(z, 'survfit'))))
            stop("input must be a square matrix survival curves")
    if (sum(nonzero) < 2)
        stop("input must have at least 2 transitions")
 
    classes <- lapply(curves, class)
    # Make sure we were sent the right things.  If any of the curves inherit
    #  from "survfitms", they are the result of a prior AJ computation;
    #  such recursion does not lead to valid estimates.
    dd <- dim(curves[[1]])
    ncurve <- prod(dd)
    for (i in 1:length(classes)) {
        if (length(classes[[i]]) != length(classes[[1]]) ||
            any(classes[[i]] != classes[[1]]))
            stop("all curves must be the same type")
        if (length(dim(curves[[i]])) != length(dd) ||
            any(dim(curves[[i]]) != dd))
            stop("all curves must be of the same dimension")
    }
    if (any(sapply(curves, function(x) inherits(x, "survfitms"))))
        stop("multi-state curves are not a valid input")
    type <- classes[[1]][1]  # 'survfit' or 'survfit.cox'

    temp <- sapply(curves, function(x) x$start.time)
    tlen <- sapply(temp, length)
    if (any(tlen >0)) {  # at least one curve with start.time
        if (any(tlen != 1) ||
            any(temp != temp[1]))
            stop("all curves must have a consistent start.time value")
        if (!missing(start.time) && temp[1] > start.time) 
            warning("curves have a larger start.time than the parameter, start.time parameter value was ignored")
        start.time <- temp[1]
    } 
        

    if (missing(method)) {
        if (type=='survfit.cox') method <- "matexp"
        else method <- "discrete"
    }
    else method <- match.arg(method)
    
    if (missing(p0)) p0 <- c(1, rep(0, nstate-1))
    if (!is.matrix(p0)) p0 <- matrix(rep(p0, ncurve), ncol=nstate, byrow=TRUE)
    else if (nrow(p0) != ncurve) stop("wrong number of rows for p0")
    if (ncol(p0) != nstate) stop("incorrect number of states in p0")
    if (any( rowSums(p0) !=1) || any(p0<0))
        stop("invalid elements in p0")

    docurve <- function(z, nonzero, nstate, p0) {
        # z is a list of survival curves
        utime <- lapply(z, function(x) x$time[x$n.event>0])
        utime <- sort(unique(unlist(utime)))  # set of unique times
        cumhaz<- lapply(z, function(x)
                    summary(x, times=utime, extend=TRUE)$cumhaz)
        jumps <- matrix(unlist(lapply(cumhaz, function(x) diff(c(0, x)))),
                        ncol= sum(nonzero))
        Tmat <- diag(nstate)
        
        # deal with the start time argument
        if (is.null(start.time)) stime <- min(c(0, utime))
        else stime <- start.time
        toss <- (utime < stime) 
        if (any(toss)) {
            utime <- utime[!toss]
            jumps <- jumps[!toss,]
            } 
        if (utime[1] > stime) {
            utime <- c(stime, utime)
            i2 <- 2
        }  else i2 <- 1
        pstate  <- matrix(0., nrow= length(utime), ncol=nstate)
        pstate[1,] <- ptemp <- p0
        for (i in i2:length(utime)) {
            Tmat[nonzero] <- jumps[i,]
            if (method == "matrix") {
                temp <- pmin(1, rowSums(Tmat) - diag(Tmat)) # failsafe
                diag(Tmat) <- 1 - temp  #rows sum to 1
                ptemp <- ptemp %*% Tmat
            }
            else {
                diag(Tmat) <- diag(Tmat) - rowSums(Tmat) #rows sum to 0
                ptemp <- as.vector(ptemp %*% expm(Tmat))
            }
            pstate[i,] <- ptemp
        }

        # Fill in the n.risk and n.event matrices
        zz <- matrix(0, nstate, nstate)
        from <- (row(zz))[nonzero]
        to   <- (col(zz))[nonzero]
        n.risk <- n.event <- matrix(0, length(utime), ncol= nstate)
        # the n.risk matrix is based on "from", n.event on "to"
        # If multiple curves come from the same source, we blithely
        #  assume that they will agree on the sample size.  If multiples
        #  go to the same ending, add the events.
        for (i in 1:length(z)) {
            index <- findInterval(utime, z[[i]]$time)
            n.risk[,from[i]] <- c(0, z[[i]]$n.risk)[index +1]
            n.event[, to[i]] <- n.event[,to[i]] + c(0, z[[i]]$n.event)[index+1]
        }
        # All the curves should have the same n
        list(n = z[[1]]$n, time = utime, pstate= pstate[-1,], 
             n.risk= n.risk, n.event=n.event)
    }
        
    # The output will have nstate columns, one for each state, and 
    #  prod(dim(curves)) strata.  If the input has strata and
    #  columns, the output strata will repeat the original strata, once
    #  for each column of the input, adding "new1", "new2", etc to the
    #  end to recognize the rows of the newdata frame that gave rise to it.
    #
    nstrat <- length(curves[[1]]$strata)
    tlist <- vector("list", prod(dd))  
    k <- 1
    if (length(dd) ==1) {
        tname <- paste0("new", 1:dd[1])
        for (j in 1:dd[1]) {
            tlist[[j]] <- docurve(lapply(curves, function(x) x[j]),
                                  nonzero, nstate, p0[j,])
        }
    } else {
        tname <-  paste0("new", 1:dd[2])
        for (j in 1:dd[2]) {
            for (i in 1:dd[1]) {
                tlist[[k]] <- docurve(lapply(curves, function(x) x[i,j]),
                                      nonzero, nstate, p0[k,])
                k <- k+1
            }
        }
    }
        
    fit <- list()
    fit$n <- tlist[[1]]$n
    fit$time <- unlist(lapply(tlist, function(x) x$time))
    fit$pstate <- do.call("rbind", lapply(tlist, function(x) x$pstate))
    fit$n.risk <- do.call("rbind", lapply(tlist, function(x) x$n.risk))
    fit$n.event<- do.call("rbind", lapply(tlist, function(x) x$n.event))
    ntemp <- unlist(lapply(tlist, function(x) length(x$time)))

    if (nstrat > 0 && length(dd)==2) 
        names(ntemp) <- as.vector(outer(names(curves[[1]]$strata), tname,
                                        paste, sep=", "))
    else names(ntemp) <- tname
    fit$strata <- ntemp
    ns <- length(fit$strata)
    fit$p0 <- p0
    fit$states <- states
    fit$n <- curves[[1]]$n
    if (!is.null(start.time)) fit$start.time <- start.time
    fit$call <- Call
    class(fit) <- c("survfitms", "survfit")
    fit
}
