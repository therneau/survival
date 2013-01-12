#
# quantile function for survfit objects
#
#  First a little function to find quantiles in a survival
#  curve.  It would be a trivial use of approx, except that
#  once in a while the survival curve has a flat spot exactly
#  at the requested quantile.  Then we use the median of the
#  flat.  
doquant <- function(p, time, surv, upper, lower, firstx) {
    findq <- function(x, y, p) {
        # special case that arises with a survival curve of only 1 subject
        if (all(is.na(y))) return(rep(NA, length(p)))
        # This one shows up with upper limits that never drop below 1.0
        if (min(y, na.rm=T) > max(1-p)) return(rep(NA, length(p)))
        
        #remove dups (censors)
        toss <- duplicated(y)
        if (any(toss)) {
            newy <- y[!toss]
            newx <- x[!toss]
        }
        else {
            newy <- y
            newx <- x
        }

        indx <- approx(1-newy, 1:length(newy), p)$y
        indx2 <- ceiling(indx)
        # quantile = where a horzontal line at p intercects the curve.  At each
        #  x the curve of 1-y jumps up to a new level
        # 
        # At this point indx != indx2 for most or all points.
        #  That is the simple case.        
        result <- c(NA, newx)[1+indx2]
        if (any(!is.na(indx) & indx==indx2)) {
            # the percentile is exactly one of the horzontal line on the graph
            # the exceptional special case is when it's the lowest one, in 
            #   which case we go all the way to the right, including censored
            #   observations
            special <- which(indx==indx2)
            upper <- c(newx, max(x))[indx2[special] +1]
            result[special] <- (result[special] + upper)/2
        }
        result
    }

    qq <- findq(c(firstx,time), c(1, surv), p) 
#   browser()
    if (missing(upper)) qq
    else rbind(qq, findq(c(firstx, time), c(1, lower), p), 
                   findq(c(firstx, time), c(1, upper), p))
    }

quantile.survfit <- function(x, probs=c(.25, .5, .75), conf.int=TRUE, ...) {
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
                                    x$lower[,i], xmin)
                    qmat[i,] <- temp[1,]
                    qupper[i,] <- temp[3,]
                    qlower[i,] <- temp[2,]
                }
                return(list(quantile=qmat, lower=qlower, upper=qupper))
            }
            else {
                for (i in 1:ncol(x$surv)) 
                    qmat[i,] <- doquant(probs, x$time, x$surv[,i], firstx=xmin)
                return(qmat)
            }
        }
        else {  
            # No strata and no matrix
            if (conf.int) {
                temp <- doquant(probs, x$time, x$surv, x$upper, x$lower, xmin)
                dimnames(temp) <- list(NULL, pname)
                return(list(quantile=temp[1,], lower=temp[2,], upper=temp[3,]))
            }
            else {
                temp <- doquant(probs, x$time, x$surv, firstx=xmin)
                names(temp) <- pname
                return(temp)
            }
        }
    }
    
    else {  
        nstrat <- length(x$strata)
        if (is.matrix(x$surv) && ncol(x$surv) >1) {
            # the least common case
            # return an array with strata as the first dimension, and
            #   the probabilites as the third.
            qmat <- array(0., dim=c(nstrat, ncol(x$surv), nprob))
            dimnames(qmat) <-list(names(x$strata), dimnames(x$surv)[[2]], pname)
            for (strat in 1:nstrat) {
                z <- x[strat]
                if (conf.int) {
                    qupper <- qlower <- qmat
                    for (i in 1:ncol(z$surv)) {
                        temp <- doquant(probs, z$time, z$surv[,i], 
                                        z$upper[,i], z$lower[,i], xmin)
                        qmat[strat,i,] <- temp[1,]
                        qupper[strat,i,] <- temp[3,]
                        qlower[strat,i,] <- temp[2,]
                    }
                }
                else {
                    for (i in 1:ncol(z$surv)) 
                        qmat[strat,i,] <- doquant(probs, z$time, z$surv[,i],
                                                  firstx=xmin)
                }
            }

            if (conf.int)                 
                    return(list(quantile=qmat, lower=qlower, upper=qupper))
            else return(qmat)
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
                                    xmin)
                    qmat[i,] <- temp[1,]
                    qupper[i,] <- temp[3,]
                    qlower[i,] <- temp[2,]
                }
                return(list(quantile=qmat, lower=qlower, upper=qupper))
            }
            else {
                for (i in 1:nstrat) {
                    z <- x[i]
                    qmat[i,] <- doquant(probs, z$time, z$surv, firstx=xmin)
                }
                return(qmat)
            }
        }
    }
}

# Why can't I just fudge the object and call quantile.survfit?  Because
#  the code below uses subscripted objects, and the class of the chimeric
#  object doesn't work out for that operation.
#
quantile.survfitms <- function(x, probs=c(.25, .5, .75), conf.int=TRUE, ...) {
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
            dimnames(qmat) <- list(dimnames(x$prev)[[2]], pname)
            if (conf.int) {
                qupper <- qlower <- qmat
                for (i in 1:ncol(x$prev)) {
                    temp <- doquant(probs, x$time, 1-x$prev[,i], 1-x$lower[,i],
                                    1-x$upper[,i], xmin)
                    qmat[i,] <- temp[1,]
                    qupper[i,] <- temp[3,]
                    qlower[i,] <- temp[2,]
                }
                return(list(quantile=qmat, lower=qlower, upper=qupper))
            }
            else {
                for (i in 1:ncol(x$prev)) 
                    qmat[i,] <- doquant(probs, x$time, 1-x$prev[,i], 
                                        firstx=xmin)
                return(qmat)
            }
        }
        else {  
            # No strata and no matrix
            if (conf.int) {
                temp <- doquant(probs, x$time, 1-x$prev,  1-x$lower, 1-x$upper,
                                firstx=xmin)
                dimnames(temp) <- list(NULL, pname)
                return(list(quantile=temp[1,], lower=temp[2,], upper=temp[3,]))
            }
            else {
                temp <- doquant(probs, x$time, 1-x$prev, firstx=xmin)
                names(temp) <- pname
                return(temp)
            }
        }
    }
    
    else {  
        nstrat <- length(x$strata)
        if (is.matrix(x$prev) && ncol(x$prev) >1) {
            # the least common case
            # return an array with strata as the first dimension, and
            #   the probabilites as the third.
            qmat <- array(0., dim=c(nstrat, ncol(x$prev), nprob))
            dimnames(qmat) <-list(names(x$strata), dimnames(x$prev)[[2]], pname)
            for (strat in 1:nstrat) {
                z <- x[strat]
                if (conf.int) {
                    qupper <- qlower <- qmat
                    for (i in 1:ncol(z$prev)) {
                        temp <- doquant(probs, z$time, 1-z$prev[,i], 
                                        1-z$lower[,i], 1-z$upper[,i], xmin)
                        qmat[strat,i,] <- temp[1,]
                        qupper[strat,i,] <- temp[3,]
                        qlower[strat,i,] <- temp[2,]
                    }
                }
                else {
                    for (i in 1:ncol(z$prev)) 
                        qmat[strat,i,] <- doquant(probs, z$time, 1-z$prev[,i],
                                                  firstx= xmin)
                }
            }

            if (conf.int)                 
                    return(list(quantile=qmat, lower=qlower, upper=qupper))
            else return(qmat)
        }       
        else {  
            # Only a strata, the most common case
            qmat <- matrix(0., nstrat, nprob)
            dimnames(qmat) <- list(names(x$strata), pname)
            if (conf.int) {
                qupper <- qlower <- qmat
                for (i in 1:nstrat) {
                    z <- x[i]
                    temp <- doquant(probs, z$time, 1-z$prev, 1-z$lower,
                                    1-z$upper, xmin)
                    qmat[i,] <- temp[1,]
                    qupper[i,] <- temp[3,]
                    qlower[i,] <- temp[2,]
                }
                return(list(quantile=qmat, lower=qlower, upper=qupper))
            }
            else {
                for (i in 1:nstrat) {
                    z <- x[i]
                    qmat[i,] <- doquant(probs, z$time, 1-z$prev, firstx=xmin)
                }
                return(qmat)
            }
        }
    }
}
