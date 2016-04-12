# Create the Aalen-Johansen estimate by joining a set of 
#  survival curves.  Actually the cumulative hazard estimates are used.
#
survfit.matrix <- function(formula, p0, method=c("discrete", "matexp"), ...) {
    Call <- match.call()
    curves <- formula
    if (!is.matrix(curves)) 
             stop("input must be a square matrix of survival curves")
    if (!is.list(curves)) stop("input must be a matrix of survival curves")
    if (nrow(curves) != ncol(curves)) 
            stop("input must be a square matrix survival curves")
    nstate <- nrow(curves)

    if (missing(p0)) p0 <- c(1, rep(0, nstate-1))
    else {
        if (sum(p0) !=1 || any(p0 <0) || length(p0) != nstate)
            stop("invalid p0 vector")
        if (length(p0) != nstate) 
            stop("input matrix and p0 have different dimensions")
    }

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

    if (missing(method)) {
        if (type=='survfit.cox') method <- "matexp"
        else method <- "discrete"
    }
    else method <- match.arg(method)
    
    docurve <- function(z, nonzero, nstate) {
        # z is a list of survival curves
        utime <- lapply(z, function(x) x$time[x$n.event>0])
        utime <- sort(unique(unlist(utime)))  # set of unique times
        cumhaz<- lapply(z, function(x)
                    summary(x, times=utime, extend=TRUE)$cumhaz)
        jumps <- matrix(unlist(lapply(cumhaz, function(x) diff(c(0, x)))),
                        ncol= sum(nonzero))
        Tmat <- diag(nstate)
        prev  <- matrix(0., nrow= 1+length(utime), ncol=nstate)
        prev[1,] <- p0
        for (i in 1:length(utime)) {
            Tmat[nonzero] <- jumps[i,]
            if (method == "matrix") {
                temp <- pmin(1, rowSums(Tmat) - diag(Tmat)) # failsafe
                diag(Tmat) <- 1 - temp  #rows sum to 1
                prev[i+1,] <- prev[i,] %*% Tmat
            }
            else {
                diag(Tmat) <- diag(Tmat) - rowSums(Tmat) #rows sum to 0
                prev[i+1,] <- as.vector(prev[i,] %*% expm(Tmat))
            }
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
            n.event[, to[i]] <- n.event[,to[i]] + z[[i]]$n.event[index]
        }
        # All the curves should have the same n
        list(n = z[[1]]$n, time = utime, prev= prev[-1,], 
             n.risk= n.risk, n.event=n.event)
    }
        
    # The output will have nstate columns, one for each state, and 
    #  prod(dim(curves)) strata.  If the input has strata and
    #  columns, the output strata will repeat the original strata, once
    #  for each column of the input, adding "row1", "row2", etc to the
    #  end to recognize the rows of the newdata frame that gave rise to it.
    #
    nstrat <- length(curves[[1]]$strata)
    if (length(dd) ==1) {
        if (dd==1) {  #only one curve
            fit <- docurve(curves, nonzero, nstate)
        } else {
            # the most common case is multiple target values in newdata
            #  but no strata
            tlist <- vector("list", dd)
            for (i in 1:dd)
                tlist[[i]] <- docurve(lapply(curves, function(x) x[i]),
                                      nonzero, nstate)
            fit <- list()
            fit$n <- tlist[[1]]$n
            fit$time <- unlist(lapply(tlist, function(x) x$time))
            fit$prev <- do.call("rbind", lapply(tlist, function(x) x$prev))
            fit$n.risk <- do.call("rbind", lapply(tlist, function(x) x$n.risk))
            fit$n.event<- do.call("rbind", lapply(tlist, function(x) x$n.event))
            ntemp <- unlist(lapply(tlist, function(x) length(x$time)))
            if (nstrat > 0) names(ntemp) <- names(strata)
            else  names(ntemp) <- paste0("new", 1:dd)
            fit$strata <- ntemp
        }
    }
    else {  
        # both strata and multiple columns
        tlist <- vector("list", prod(dd))
        k <- 1
        for (j in 1:dd[2]) {
            for (j in 1:dd[1]) {
                tlist[[k]] <- docurve(lapply(curves, function(x) x[i,j]),
                                      nonzero, nstate)
                k <- k+1
            }
        }
        
        fit <- list()
        fit$n <- tlist[[1]]$n
        fit$time <- unlist(lapply(tlist, function(x) x$time))
        fit$prev <- do.call("rbind", lapply(tlist, function(x) x$prev))
        fit$n.risk <- do.call("rbind", lapply(tlist, function(x) x$n.risk))
        fit$n.event<- do.call("rbind", lapply(tlist, function(x) x$n.event))
        ntemp <- unlist(lapply(tlist, function(x) length(x$time)))
        names(ntemp) <- as.vector(outer(names(strata), paste0("new", 1:dd[2]), 
                                        paste, sep=", "))
        fit$strata <- ntemp
    }

    if (length(fit$strata) > 0) {
        ns <- length(fit$strata)
        fit$p0 <- matrix(rep(p0, each=ns), nrow=ns)
    }
    else fit$p0 <- p0
    fit$states <- states
    fit$n <- curves[[1]]$n
    fit$call <- Call
    class(fit) <- c("survfitms", "survfit")
    fit
}
