#
# quantile function for survfit objects
#
#  First a little function to find quantiles in a CDF
#  curve.  It would be a trivial use of approx, except that
#  once in a while the survival curve has a flat spot exactly
#  at the requested quantile.  Then we use the median of the
#  flat.  
findq <- function(x, y, p, tol) {
    # This case occurs for a survival curve whose upper limit never drops below 1
    if (max(y, na.rm=T) < min(p)) return(rep(NA, length(p)))
        
    # Remove duplicate y values, i.e., the censors, since dups cause
    #  issues for approx
    xmax <- x[length(x)]
    dups <- duplicated(y)
    if (any(dups)) {
        x <- x[!dups]
        y <- y[!dups]
    }
    n <- length(y)

    # quantile = where a horzontal line at p intercects the curve.  At each
    #  x the curve of 1-y jumps up to a new level
    # The most work is to check for horizontal lines in the survival 
    #  curve that match one of our quantiles within tolerance.  If any
    #  p matches, then our quantile is the average of the given x and
    #  the x value of the next jump point, i.e., the usual midpoint rule
    #  used for medians.  
    # A flat at the end of the curve is a special case, as is the quantile
    #  of 0.
    indx1 <- approx(y+tol, 1:n, p, method="constant", f=1)$y
    indx2 <- approx(y-tol, 1:n, p, method="constant", f=1)$y
    quant <- (x[indx1] + x[indx2])/2
    quant[p==0] <- x[1]
    if (!is.na(y[n])) {
        lastpt <- (abs(p- y[n]) < tol)  # end of the curve
        if (any(lastpt)) quant[lastpt] <- (x[indx1[lastpt]] + xmax)/2
    }
    quant
    }

doquant <- function(p, time, surv, upper, lower, firstx, tol) {
    qq <- findq(c(firstx,time), c(0, 1-surv), p, tol) 
#   browser()
    if (missing(upper)) qq
    else rbind(qq, findq(c(firstx, time), c(0, 1-lower), p, tol), 
                   findq(c(firstx, time), c(0, 1-upper), p, tol))
    }

quantile.survfit <- function(x, probs=c(.25, .5, .75), conf.int=TRUE, 
                             tolerance= sqrt(.Machine$double.eps), ...) {
    if (!inherits(x, "survfit")) stop("Must be a survfit object")
    if (any(!is.numeric(probs)) || any(is.na(probs)))
        stop("invalid probability")
    if (any(probs <0 | probs >1)) stop("Invalid probability")
    if (is.null(x$lower)) conf.int <- FALSE
    nprob <- length(probs)
    pname <- format(probs*100)

    # What do we report for p=0?  Use x$start.time if it exists, 0 otherwise
    xmin <- if (is.null(x$start.time)) 0 else x$start.time

    # There are 8 cases: strata yes/no
    #                    ncol(x$surv) =1 or >1
    #                    conf.int = T/F
    if (is.null(x$strata)) {
        if (is.matrix(x$surv) && ncol(x$surv) >1) {
            qmat <- matrix(0., ncol=nprob, nrow=ncol(x$surv))
            dimnames(qmat) <- list(dimnames(x$surv)[[2]], pname)
            if (conf.int) {
                qupper <- qlower <- qmat
                for (i in 1:ncol(x$surv)) {
                    temp <- doquant(probs, x$time, x$surv[,i], x$upper[,i],
                                    x$lower[,i], xmin, tolerance)
                    qmat[i,] <- temp[1,]
                    qupper[i,] <- temp[3,]
                    qlower[i,] <- temp[2,]
                }
                return(list(quantile=qmat, lower=qlower, upper=qupper))
            }
            else {
                for (i in 1:ncol(x$surv)) 
                    qmat[i,] <- doquant(probs, x$time, x$surv[,i], firstx=xmin,
                                        tol=tolerance)
                return(qmat)
            }
        }
        else {  
            # No strata and no matrix
            if (conf.int) {
                temp <- doquant(probs, x$time, x$surv, x$upper, x$lower, xmin,
                                tolerance)
                dimnames(temp) <- list(NULL, pname)
                return(list(quantile=temp[1,], lower=temp[2,], upper=temp[3,]))
            }
            else {
                temp <- doquant(probs, x$time, x$surv, firstx=xmin, 
                                tol =tolerance)
                names(temp) <- pname
                return(temp)
            }
        }
    }
    
    else {  
        nstrat <- length(x$strata)
        if (is.matrix(x$surv) && ncol(x$surv) >1) {
            # uncommon case, e.g., predicted survivals from a Cox model
            # return an array with strata as the first dimension, and
            #   the probabilites as the third.
            qmat <- array(0., dim=c(nstrat, ncol(x$surv), nprob))
            dimnames(qmat) <-list(names(x$strata), dimnames(x$surv)[[2]], pname)
            if (conf.int) {
                qupper <- qlower <- qmat
                for (strat in 1:nstrat) {
                    z <- x[strat,]
                    for (i in 1:ncol(z$surv)) {
                        temp <- doquant(probs, z$time, z$surv[,i], 
                                        z$upper[,i], z$lower[,i], xmin,tolerance)
                        qmat[strat,i,] <- temp[1,]
                        qupper[strat,i,] <- temp[3,]
                        qlower[strat,i,] <- temp[2,]
                    }
                }
                return(list(quantile=qmat, lower=qlower, upper=qupper))
            }
            else {
                for (strat in 1:nstrat) {
                    z <- x[strat]
                    for (i in 1:ncol(z$surv)) 
                        qmat[strat,i,] <- doquant(probs, z$time, z$surv[,i],
                                                  firstx=xmin, tol=tolerance)
                }
                return(qmat)
            }
         }       
        else {  
            # Only a strata, the most common case
            qmat <- matrix(0., nstrat, nprob)
            dimnames(qmat) <- list(names(x$strata), pname)
            if (conf.int) {
                qupper <- qlower <- qmat
                for (i in 1:nstrat) {
                    z <- x[i]
                    temp <- doquant(probs, z$time, z$surv, z$upper, z$lower,
                                    xmin, tolerance)
                    qmat[i,] <- temp[1,]
                    qupper[i,] <- temp[3,]
                    qlower[i,] <- temp[2,]
                }
                return(list(quantile=qmat, lower=qlower, upper=qupper))
            }
            else {
                for (i in 1:nstrat) {
                    z <- x[i]
                    qmat[i,] <- doquant(probs, z$time, z$surv, firstx=xmin,
                                        tol = tolerance)
                }
                return(qmat)
            }
        }
    }
}

# Why can't I just fudge the object and call quantile.survfit?  Because
#  the code below uses subscripted objects, and the class of the chimeric
#  object doesn't work out for that operation.  Also, we want to use the
#  state names in the dimnames of the result.
#
quantile.survfitms <- function(x, probs=c(.25, .5, .75), conf.int=TRUE, 
                               tolerance= sqrt(.Machine$double.eps), ...) {
    if (any(!is.numeric(probs)) || any(is.na(probs)))
        stop("invalid probability")
    if (any(probs <0 | probs >1)) stop("Invalid probability")
    if (is.null(x$lower)) conf.int <- FALSE
    nprob <- length(probs)
    pname <- format(probs*100)

    if (is.null(x$start.time)) xmin<-0 else xmin <- x$start.time

    # There are 8 cases: strata yes/no
    #                    ncol(x$surv) =1 or >1
    #                    conf.int = T/F
    if (is.null(x$strata)) {
        if (is.matrix(x$prev) && ncol(x$prev) >1) {
            qmat <- matrix(0., ncol=nprob, nrow=ncol(x$prev))
            dimnames(qmat) <- list(x$states, pname)
            if (conf.int) {
                qupper <- qlower <- qmat
                for (i in 1:ncol(x$prev)) {
                    temp <- doquant(probs, x$time, 1-x$prev[,i], 1-x$lower[,i],
                                    1-x$upper[,i], xmin, tolerance)
                    qmat[i,] <- temp[1,]
                    qupper[i,] <- temp[3,]
                    qlower[i,] <- temp[2,]
                }
                return(list(quantile=qmat, lower=qlower, upper=qupper))
            }
            else {
                for (i in 1:ncol(x$prev)) 
                    qmat[i,] <- doquant(probs, x$time, 1-x$prev[,i], 
                                        firstx=xmin, tol=tolerance)
                return(qmat)
            }
        }
        else {  
            # No strata and no matrix
            if (conf.int) {
                temp <- doquant(probs, x$time, 1-x$prev,  1-x$lower, 1-x$upper,
                                firstx=xmin, tol=tolerance)
                dimnames(temp) <- list(NULL, pname)
                return(list(quantile=temp[1,], lower=temp[2,], upper=temp[3,]))
            }
            else {
                temp <- doquant(probs, x$time, 1-x$prev, firstx=xmin, 
                                tol=tolerance)
                names(temp) <- pname
                return(temp)
            }
        }
    }
    
    else {  
        nstrat <- length(x$strata)
        if (is.matrix(x$prev) && ncol(x$prev) >1) {
            # the most common case
            # Return an array with strata as the first dimension, state as
            #   the second, and the probabilites as the third.  The reason
            #   for this order is that then 
            #      (quantile(fit))[i,j,] = quantile(fit[i,j])
            qmat <- array(0., dim=c(nstrat, ncol(x$prev), nprob))
            dimnames(qmat) <-list(names(x$strata), x$states, pname)
            if (conf.int) {
                qupper <- qlower <- qmat
                for (strat in 1:nstrat) {
                    z <- x[strat]
                    for (i in 1:ncol(z$prev)) {
                        temp <- doquant(probs, z$time, 1-z$prev[,i], 
                                 1-z$lower[,i], 1-z$upper[,i], xmin, tolerance)
                        qmat[strat,i,] <- temp[1,]
                        qupper[strat,i,] <- temp[3,]
                        qlower[strat,i,] <- temp[2,]
                    }
                }
                return(list(quantile=qmat, lower=qlower, upper=qupper))
            }
            else {
                for (strat in 1:nstrat) {
                    z <- x[strat]
                    for (i in 1:ncol(z$prev)) 
                        qmat[strat,i,] <- doquant(probs, z$time, 1-z$prev[,i],
                                                  firstx= xmin, tol=tolerance)
                }
                return(qmat)
            }
        }       
        else {  
            # Only a strata, which will be a rare case.  Perhaps someone
            #  typed "quantile(fit[,2], .4)"  
            qmat <- matrix(0., nstrat, nprob)
            dimnames(qmat) <- list(names(x$strata), pname)
            if (conf.int) {
                qupper <- qlower <- qmat
                for (i in 1:nstrat) {
                    z <- x[i]
                    temp <- doquant(probs, z$time, 1-z$prev, 1-z$lower,
                                    1-z$upper, xmin, tolerance)
                    qmat[i,] <- temp[1,]
                    qupper[i,] <- temp[3,]
                    qlower[i,] <- temp[2,]
                }
                return(list(quantile=qmat, lower=qlower, upper=qupper))
            }
            else {
                for (i in 1:nstrat) {
                    z <- x[i]
                    qmat[i,] <- doquant(probs, z$time, 1-z$prev, firstx=xmin,
                                        tol = tolerance)
                }
                return(qmat)
            }
        }
    }
}
