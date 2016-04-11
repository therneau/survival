# Create the Aalen-Johansen estimate by joining a set of 
#  survival curves.  Actually the cumulative hazard estimates are used.
#
survfit.matrix <- function(x, p0, method=c("discrete", "matexp")) {
    Call <- match.call()
    curves <- x
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
    
    docurve <- function(z) {
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
                prev[i+1,] <- prev[i,] %*% expm(Tmat)
            }
        }
        list(time = utime, pstate= prev[-1,])
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
            fit <- docurve(curves)
            fit$n.risk <- curves[[1]]$n.risk
            fit$n.event<- curves[[1]]$n.event
        }
        if (nstrat ==0) { 
            # the most common case: multiple target values in newdata
            #  but no strata
            temp <- vector("list", dd)
            for (i in 1:dd)
                temp[[i]] <- docurve(lapply(curves, function(x) x[i]))
            nn <- sapply(temp, function(x) length(x$time))
            if (any(nn != nn[1])) stop("internal error 1, survfit.matrix")
            fit <- list()
            fit$time <- unlist(sapply(temp, function(x) x$time))
            fit$pstate <- do.call("rbind", lapply(temp, function(x) x$pstate))
            index <- match(fit$time, curves[[1]]$time)
            fit$n.risk <- curves[[1]]$n.risk[index]
            fit$n.event<- curves[[1]]$n.event[index]
            names(nn) <- paste0("row", 1:dd)
            fit$strata <- nn
        } else {
            # one target value in newdata, multiple strata in the Cox model
            temp <- vector("list", dd)
            for (i in 1:dd)
                temp[[i]] <- docurve(lapply(curves, function(x) x[i]))
            nn <- sapply(temp, function(x) length(x$time))
            fit <- list()
            fit$time <- unlist(sapply(temp, function(x) x$time))
            fit$pstate <- do.call("rbind", sapply(temp, function(x) x$pstate))
            fit$n.risk <- unlist(sapply(1:dd, function(i) {
                index <- match(temp[[i]]$time, curves[[1]][i]$time)
                curves[[1]][i]$n.risk[index]
            }))
            fit$n.event <- unlist(sapply(1:dd, function(i) {
                index <- match(temp[[i]]$time, curves[[1]][i]$time)
                curves[[1]][i]$n.event[index]
            }))
            names(nn) <- names(curves[[1]]$strata)
            fit$strata <- nn
        }
    }
    else{  # both strata and newdata, yipes!
        stop("code not yet finished for this case ")
    }
    if (length(fit$strata) > 0) {
        ns <- length(fit$strata)
        fit$p0 <- matrix(rep(p0, each=ns), nrow=ns)
    }
    else fit$p0 <- p0
    fit$states <- states
    fit$n <- curves[[1]]$n
    class(fit) <- c("survfitms", "survfit")
    fit
}
