#
# quantile function for survfit objects
#
#  First a little function to find quantiles in a survival
#  curve.  It would be a trivial use of approx, except that
#  one in a while the survival curve has a flat spot exactly
#  at the requested quantile.  Then we use the median of the
#  flat.  
doquant <- function(p, time, surv, upper, lower) {
    findq <- function(x, y, p) {
        # special case that arises with a survival curve of only 1 subject
        if (all(is.na(y))) return(rep(NA, length(p)))
        
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

        indx <- approx(c(0, 1-newy), c(0:length(newy)), p)$y 
        indx2 <- ceiling(indx)
        # quantile = where a horzontal line at p intercects the curve.  At each
        #  x the curve of 1-y jumps up to a new level
        # 
        # At this point indx != indx2 for most or all points.
        #  That is the simple case.        
        result <- newx[indx2]
        if (any(!is.na(indx) & indx==indx2)) {
            # the percentile is exactly one of the horzontal line on the graph
            # the exceptional special case is when it's the lowest one, in 
            #   which case we go all the way to the right, including censored
            #   observations
            special <- which(indx==indx2)
            upper <- c(new, max(x))[indx2[special] +1]
            result[special] <- (result[special] + upper)/2
        }
        result
    }

    qq <- findq(time, surv, p)
    if (missing(upper)) qq
    else rbind(qq, findq(time, lower, p), findq(time, upper, p))
    }

quantile.survfit <- function(x, probs=c(.25, .5, .75), conf.int=TRUE, ...) {
    if (!inherits(x, "survfit")) stop("Must be a survfit object")

    if (any(probs <0 | probs >=1)) stop("Invalid probability")
    if (is.null(x$lower)) conf.int <- FALSE
    nprob <- length(probs)
    pname <- paste(probs*100, "%", sep="")

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
                                    x$lower[,i])
                    qmat[i,] <- temp[1,]
                    qupper[i,] <- temp[3,]
                    qlower[i,] <- temp[2,]
                }
                return(list(quantile=qmat, lower=qlower, upper=qupper))
            }
            else {
                for (i in 1:ncol(x$surv)) 
                    qmat[i,] <- doquant(probs, x$time, x$surv[,i])
                return(qmat)
            }
        }
        else {  
            # No strata and no matrix
            if (conf.int) {
                temp <- doquant(x$time, x$surv, x$upper, x$lower)
                dimnames(temp) <- list(NULL, pname)
                return(list(quantile=temp[1,], lower=temp[2,], upper=temp[3,]))
            }
        }
    }
    
    else {  
        nstrat <- length(x$strata)
        if (is.matrix(x$surv) && ncol(x$surv) >1) {
            # the least common case
            # return an array with strata as the 3rd dimension
            qmat <- matrix(0., ncol=nprob, nrow=ncol(x$surv))
            dimnames(qmat) <- list(dimnames(x$surv)[[2]], pname, 
                                   names(x$strata))
            for (strat in 1:nstrat) {
                z <- x$surv[strat,]
                if (conf.int) {
                    qupper <- qlower <- qmat
                    for (i in 1:ncol(z$surv)) {
                        temp <- doquant(probs, z$time, z$surv[,i], z$upper[,i],
                                    z$lower[,i])
                        qmat[i,,strat] <- temp[1,]
                        qupper[i,,strat] <- temp[3,]
                        qlower[i,,strat] <- temp[2,]
                    }
                }
                else {
                    for (i in 1:ncol(x$surv)) 
                        qmat[i,] <- doquant(probs, x$time, x$surv[,i])
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
                    temp <- doquant(probs, z$time, z$surv, z$upper, z$lower)
                    qmat[i,] <- temp[1,]
                    qupper[i,] <- temp[3,]
                    qlower[i,] <- temp[2,]
                }
                return(list(quantile=qmat, lower=qlower, upper=qupper))
            }
            else {
                for (strata in 1:nstrat) {
                    z <- x[i]
                    qmat[i,] <- doquant(z$time, z$surv)
                }
                return(qmat)
            }
        }
    }
}

