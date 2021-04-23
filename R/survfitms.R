# Automatically generated from the noweb directory
# Methods for survfitms objects
summary.survfit <- function(object, times, censored=FALSE, 
                            scale=1, extend=FALSE, 
                            rmean=getOption('survfit.rmean'),
                            ...) {
    fit <- object  #save typing
    if (!inherits(fit, 'survfit'))
            stop("summary.survfit can only be used for survfit objects")
    if (is.null(fit$logse)) fit$logse <- TRUE   #older style

    # The print.rmean option is depreciated, it is still listened
    #   to in print.survfit, but ignored here
    if (is.null(rmean)) rmean <- "common"
    if (is.numeric(rmean)) {
        if (is.null(fit$start.time)) {
            if (rmean < min(fit$time)) 
                stop("Truncation point for the mean time in state is < smallest survival")
        }
        else if (rmean < fit$start.time)
            stop("Truncation point for the mean time in state is < smallest survival")
    }
    else {
        rmean <- match.arg(rmean, c('none', 'common', 'individual'))
        if (length(rmean)==0) stop("Invalid value for rmean option")
    }

    # adding time 0 makes the mean and median easier
    fit0 <- survfit0(fit, fit$start.time)  #add time 0
    temp <- survmean(fit0, scale=scale, rmean)  
    table <- temp$matrix  #for inclusion in the output list
    rmean.endtime <- temp$end.time
    
    if (!is.null(fit$strata)) {
        nstrat <-  length(fit$strata)
    }    
    delta <- function(x, indx) {  # sums between chosen times
        if (is.logical(indx)) indx <- which(indx)
        if (!is.null(x) && length(indx) >0) {
            fx <- function(x, indx) diff(c(0, c(0, cumsum(x))[indx+1]))
            if (is.matrix(x)) {
                temp <- apply(x, 2, fx, indx=indx)
                # don't return a vector when only 1 time point is given
                if (is.matrix(temp)) temp else matrix(temp, nrow=1)
            }
            else fx(x, indx)
        }
        else NULL
    }

    if (missing(times)) {
        if (!censored) {
            index <- (rowSums(as.matrix(fit$n.event)) >0)
            for (i in c("time","n.risk", "n.event", "surv", "pstate", "std.err", 
                                "upper", "lower", "cumhaz", "std.chaz")) {
                if (!is.null(fit[[i]])) {  # not all components in all objects
                    temp <- fit[[i]]
                    if (is.matrix(temp)) temp <- temp[index,,drop=FALSE]
                    else  if (!is.array(temp)) temp <- temp[index]  #simple vector
                    else temp <- temp[index,,, drop=FALSE] # 3 way
                    fit[[i]] <- temp
                }
            }

            # The n.enter and n.censor values are accumualated
            #  both of these are simple vectors
            if (is.null(fit$strata)) {
                for (i in c("n.enter", "n.censor"))
                    if (!is.null(fit[[i]]))
                        fit[[i]] <- delta(fit[[i]], index)
            }
            else {
                sindx <- rep(1:nstrat, fit$strata)
                for (i in c("n.enter", "n.censor")) {
                    if (!is.null(fit[[i]]))
                        fit[[i]] <- unlist(sapply(1:nstrat, function(j) 
                                     delta(fit[[i]][sindx==j], index[sindx==j])))
                }
                # the "factor" is needed for the case that a strata has no
                #  events at all, and hence 0 lines of output
                fit$strata[] <- as.vector(table(factor(sindx[index], 1:nstrat))) 
            }
        }
        #if missing(times) and censored=TRUE, the fit object is ok as it is
    }
    else {
        fit <- fit0
        ssub<- function(x, indx) {  #select an object and index
            if (!is.null(x) && length(indx)>0) {
                if (is.matrix(x)) x[pmax(1,indx),,drop=FALSE]
                else if (is.array(x))  x[pmax(1,indx),,,drop=FALSE]
                else x[pmax(1, indx)]
            }
            else NULL
        }
        findrow <- function(fit, times, extend) {
            if (FALSE) {
                if (is.null(fit$start.time)) mintime <- min(fit$time, 0)
                else                         mintime <- fit$start.time
                ptimes <- times[times >= mintime]
            }  else ptimes <- times[is.finite(times)]      

            if (!extend) {
                maxtime <- max(fit$time)
                ptimes <- ptimes[ptimes <= maxtime]
            }
            ntime <- length(fit$time)
            
            index1 <- findInterval(ptimes, fit$time) 
            index2 <- 1 + findInterval(ptimes, fit$time, left.open=TRUE)
            if (length(index1) ==0)
                stop("no points selected for one or more curves, consider using the extend argument")
            # The pmax() above encodes the assumption that n.risk for any
            #  times before the first observation = n.risk at the first obs
            fit$time <- ptimes

            for (i in c("surv", "pstate", "upper", "lower", "std.err", "cumhaz",
                        "std.chaz")) {
                if (!is.null(fit[[i]])) fit[[i]] <- ssub(fit[[i]], index1)
            }
            
            if (is.matrix(fit$n.risk)) {
                # Every observation in the data has to end with a censor or event.
                #  So by definition the number at risk after the last observed time
                #  value must be 0.
                fit$n.risk <- rbind(fit$n.risk,0)[index2,,drop=FALSE]
            }
            else  fit$n.risk <- c(fit$n.risk, 0)[index2]

            for (i in c("n.event", "n.censor", "n.enter"))
                fit[[i]] <- delta(fit[[i]], index1)
            fit
        }

        # For a single component, turn it from a list into a single vector, matrix
        #  or array
        unlistsurv <- function(x, name) {
            temp <- lapply(x, function(x) x[[name]])
            if (is.vector(temp[[1]])) unlist(temp)
            else if (is.matrix(temp[[1]])) do.call("rbind", temp)
            else { 
                # the cumulative hazard is the only component that is an array
                # it's third dimension is n
                xx <- unlist(temp)
                dd <- dim(temp[[1]])
                dd[3] <- length(xx)/prod(dd[1:2])
                array(xx, dim=dd)
            }
        }

        # unlist all the components built by a set of calls to findrow
        #  and remake the strata
        unpacksurv <- function(fit, ltemp) {
            keep <- c("time", "surv", "pstate", "upper", "lower", "std.err",
                      "cumhaz", "n.risk", "n.event", "n.censor", "n.enter",
                      "std.chaz")
            for (i in keep) 
                if (!is.null(fit[[i]])) fit[[i]] <- unlistsurv(ltemp, i)
            fit$strata[] <- sapply(ltemp, function(x) length(x$time))
            fit
        }
        times <- sort(times)  #in case the user forgot
        if (is.null(fit$strata)) fit <- findrow(fit, times, extend)
        else {
            ltemp <- vector("list", nstrat)
            for (i in 1:nstrat) 
                ltemp[[i]] <- findrow(fit[i], times, extend)
            fit <- unpacksurv(fit, ltemp)
        }
    }

    # finish off the output structure
    fit$table <- table
    if (length(rmean.endtime)>0  && !any(is.na(rmean.endtime[1]))) 
            fit$rmean.endtime <- rmean.endtime

    # A survfit object may contain std(log S) or std(S), summary always std(S)
    if (!is.null(fit$std.err) && fit$logse) fit$std.err <- fit$std.err * fit$surv 
 
    # Expand the strata
    if (!is.null(fit$strata)) 
        fit$strata <- factor(rep(1:nstrat, fit$strata), 1:nstrat,
                             labels= names(fit$strata))
    if (scale != 1) {
        # fix scale in the output
        fit$time <- fit$time/scale
    }

    class(fit) <- "summary.survfit"
    fit
}
summary.survfitms <- function(object, times, censored=FALSE, 
                            scale=1, extend=FALSE, 
                            rmean= getOption("survfit.rmean"),
                            ...) {

    fit <- object  # save typing
    if (!inherits(fit, 'survfitms'))
            stop("summary.survfitms can only be used for survfitms objects")
    if (is.null(fit$logse)) fit$logse <- FALSE  # older style

    # The print.rmean option is depreciated, it is still listened
    #   to in print.survfit, but ignored here
    if (is.null(rmean)) rmean <- "common"
    if (is.numeric(rmean)) {
        if (is.null(fit$start.time)) {
            if (rmean < min(fit$time)) 
                stop("Truncation point for the mean is < smallest survival")
        }
        else if (rmean < fit$start.time)
            stop("Truncation point for the mean is < smallest survival")
    }
    else {
        rmean <- match.arg(rmean, c('none', 'common', 'individual'))
        if (length(rmean)==0) stop("Invalid value for rmean option")
    }

    fit0 <- survfit0(fit, fit$start.time) # add time 0
    temp <- survmean2(fit0, scale=scale, rmean)  
    table <- temp$matrix  #for inclusion in the output list
    rmean.endtime <- temp$end.time

    if (!missing(times)) {
        if (!is.numeric(times)) stop ("times must be numeric")
        times <- sort(times)
    }

    if (!is.null(fit$strata)) {
        nstrat <-  length(fit$strata)
        sindx <- rep(1:nstrat, fit$strata)
    }    
    delta <- function(x, indx) {  # sums between chosen times
        if (is.logical(indx)) indx <- which(indx)
        if (!is.null(x) && length(indx) >0) {
            fx <- function(x, indx) diff(c(0, c(0, cumsum(x))[indx+1]))
            if (is.matrix(x)) {
                temp <- apply(x, 2, fx, indx=indx)
                if (is.matrix(temp)) temp else matrix(temp, nrow=1)
            }
            else fx(x, indx)
        }
        else NULL
    }

    if (missing(times)) {
        if (!censored) {
            index <- (rowSums(as.matrix(fit$n.event)) >0)
            for (i in c("time","n.risk", "n.event", "surv", "pstate", "std.err", 
                                "upper", "lower", "cumhaz", "std.chaz")) {
                if (!is.null(fit[[i]])) {  # not all components in all objects
                    temp <- fit[[i]]
                    if (is.matrix(temp)) temp <- temp[index,,drop=FALSE]
                    else  if (!is.array(temp)) temp <- temp[index]  #simple vector
                    else temp <- temp[index,,, drop=FALSE] # 3 way
                    fit[[i]] <- temp
                }
            }

            # The n.enter and n.censor values are accumualated
            #  both of these are simple vectors
            if (is.null(fit$strata)) {
                for (i in c("n.enter", "n.censor"))
                    if (!is.null(fit[[i]]))
                        fit[[i]] <- delta(fit[[i]], index)
            }
            else {
                sindx <- rep(1:nstrat, fit$strata)
                for (i in c("n.enter", "n.censor")) {
                    if (!is.null(fit[[i]]))
                        fit[[i]] <- unlist(sapply(1:nstrat, function(j) 
                                     delta(fit[[i]][sindx==j], index[sindx==j])))
                }
                # the "factor" is needed for the case that a strata has no
                #  events at all, and hence 0 lines of output
                fit$strata[] <- as.vector(table(factor(sindx[index], 1:nstrat))) 
            }
        }
        #if missing(times) and censored=TRUE, the fit object is ok as it is
    }
    else {
        fit <-fit0  # easier to work with
        ssub<- function(x, indx) {  #select an object and index
            if (!is.null(x) && length(indx)>0) {
                if (is.matrix(x)) x[pmax(1,indx),,drop=FALSE]
                else if (is.array(x))  x[pmax(1,indx),,,drop=FALSE]
                else x[pmax(1, indx)]
            }
            else NULL
        }
        findrow <- function(fit, times, extend) {
            if (FALSE) {
                if (is.null(fit$start.time)) mintime <- min(fit$time, 0)
                else                         mintime <- fit$start.time
                ptimes <- times[times >= mintime]
            }  else ptimes <- times[is.finite(times)]      

            if (!extend) {
                maxtime <- max(fit$time)
                ptimes <- ptimes[ptimes <= maxtime]
            }
            ntime <- length(fit$time)
            
            index1 <- findInterval(ptimes, fit$time) 
            index2 <- 1 + findInterval(ptimes, fit$time, left.open=TRUE)
            if (length(index1) ==0)
                stop("no points selected for one or more curves, consider using the extend argument")
            # The pmax() above encodes the assumption that n.risk for any
            #  times before the first observation = n.risk at the first obs
            fit$time <- ptimes

            for (i in c("surv", "pstate", "upper", "lower", "std.err", "cumhaz",
                        "std.chaz")) {
                if (!is.null(fit[[i]])) fit[[i]] <- ssub(fit[[i]], index1)
            }
            
            if (is.matrix(fit$n.risk)) {
                # Every observation in the data has to end with a censor or event.
                #  So by definition the number at risk after the last observed time
                #  value must be 0.
                fit$n.risk <- rbind(fit$n.risk,0)[index2,,drop=FALSE]
            }
            else  fit$n.risk <- c(fit$n.risk, 0)[index2]

            for (i in c("n.event", "n.censor", "n.enter"))
                fit[[i]] <- delta(fit[[i]], index1)
            fit
        }

        # For a single component, turn it from a list into a single vector, matrix
        #  or array
        unlistsurv <- function(x, name) {
            temp <- lapply(x, function(x) x[[name]])
            if (is.vector(temp[[1]])) unlist(temp)
            else if (is.matrix(temp[[1]])) do.call("rbind", temp)
            else { 
                # the cumulative hazard is the only component that is an array
                # it's third dimension is n
                xx <- unlist(temp)
                dd <- dim(temp[[1]])
                dd[3] <- length(xx)/prod(dd[1:2])
                array(xx, dim=dd)
            }
        }

        # unlist all the components built by a set of calls to findrow
        #  and remake the strata
        unpacksurv <- function(fit, ltemp) {
            keep <- c("time", "surv", "pstate", "upper", "lower", "std.err",
                      "cumhaz", "n.risk", "n.event", "n.censor", "n.enter",
                      "std.chaz")
            for (i in keep) 
                if (!is.null(fit[[i]])) fit[[i]] <- unlistsurv(ltemp, i)
            fit$strata[] <- sapply(ltemp, function(x) length(x$time))
            fit
        }
        times <- sort(times)
        if (is.null(fit$strata)) fit <- findrow(fit, times, extend)
        else {
            ltemp <- vector("list", nstrat)
            for (i in 1:nstrat) 
                ltemp[[i]] <- findrow(fit[i,], times, extend)
            fit <- unpacksurv(fit, ltemp)
        }
    }

    # finish off the output structure
    fit$table <- table
    if (length(rmean.endtime)>0  && !any(is.na(rmean.endtime))) 
            fit$rmean.endtime <- rmean.endtime

    if (!is.null(fit$strata)) 
        fit$strata <- factor(rep(names(fit$strata), fit$strata))

    # A survfit object may contain std(log S) or std(S), summary always std(S)
    if (!is.null(fit$std.err) && fit$logse) fit$std.err <- fit$std.err * fit$surv 
    if (scale != 1) {
        # fix scale in the output
        fit$time <- fit$time/scale
    }
    class(fit) <- "summary.survfitms"
    fit
}

print.survfitms <- function(x, scale=1,
                            rmean = getOption("survfit.rmean"), ...) {
    if (!is.null(cl<- x$call)) {
        cat("Call: ")
        dput(cl)
        cat("\n")
        }        
    omit <- x$na.action
    if (length(omit)) cat("  ", naprint(omit), "\n")

    x <- survfit0(x, x$start.time)
    if (is.null(rmean)) rmean <- "common"
    if (is.numeric(rmean)) {
        if (is.null(x$start.time)) {
            if (rmean < min(x$time)) 
                stop("Truncation point for the mean is < smallest survival")
        }
        else if (rmean < x$start.time)
            stop("Truncation point for the mean is < smallest survival")
    }
    else {
        rmean <- match.arg(rmean, c('none', 'common', 'individual'))
        if (length(rmean)==0) stop("Invalid value for rmean option")
    }

    temp <- survmean2(x, scale=scale, rmean)
    if (is.null(temp$end.time)) print(temp$matrix, ...)
    else {
        etime <- temp$end.time
        dd <- dimnames(temp$matrix)
        cname <- dd[[2]]
        cname[length(cname)] <- paste0(cname[length(cname)], '*')
        dd[[2]] <- cname
        dimnames(temp$matrix) <- dd
        print(temp$matrix, ...)
        if (length(etime) ==1)
             cat("   *mean time in state, restricted (max time =", 
                 format(etime, ...), ")\n")
        else cat("   *mean time in state, restricted (per curve cutoff)\n")
    }
    invisible(x)
}
survmean2 <- function(x, scale=1, rmean) {
    nstate <- length(x$states)  #there will always be at least 1 state
    ngrp   <- max(1, length(x$strata))
    if (is.null(x$newdata)) ndata <- 0  else ndata <- nrow(x$newdata)
    if (ngrp >1)  {
        igrp <- rep(1:ngrp, x$strata)
        rname <- names(x$strata)
        }
    else {
        igrp <- rep(1, length(x$time))
        rname <- NULL
        }

    # The n.event matrix may not have nstate columms.  Its
    #  colnames are the first elements of states, however
    if (is.matrix(x$n.event)) {
        nc <- ncol(x$n.event)
        nevent <- tapply(x$n.event, list(rep(igrp, nc), col(x$n.event)), sum)
        dimnames(nevent) <- list(rname, x$states[1:nc])
        }
    else {
        nevent <- tapply(x$n.event, igrp, sum)
        names(nevent) <- rname
        }

    if (ndata< 2) {
        outmat <- matrix(0., nrow=nstate*ngrp , ncol=2)
        outmat[,1] <- rep(x$n, nstate)
        outmat[1:length(nevent), 2] <- c(nevent)
        
        if (ngrp >1) 
            rowname <- c(outer(rname, x$states, paste, sep=", "))
        else rowname <- x$states
    }
    else {
        outmat <- matrix(0., nrow=nstate*ndata*ngrp, ncol=2)
        outmat[,1] <- rep(x$n, nstate*ndata)
        outmat[, 2] <- rep(c(nevent), each=ndata)
       
        temp <- outer(1:ndata, x$states, paste, sep=", ")
        if (ngrp >1) 
            rowname <- c(outer(rname, temp, paste, sep=", "))
        else rowname <- temp
        nstate <- nstate * ndata
    }

    # Caculate the mean time in each state
    if (rmean != "none") {
        if (is.numeric(rmean)) maxtime <- rep(rmean, ngrp)
        else if (rmean=="common") maxtime <- rep(max(x$time), ngrp)
        else maxtime <- tapply(x$time, igrp, max)
    
        meantime <- matrix(0., ngrp, nstate)
        if (!is.null(x$influence)) stdtime <- meantime
        for (i in 1:ngrp) {
            # a 2 dimensional matrix is an "array", but a 3-dim array is
            #  not a "matrix", so check for matrix first.
            if (is.matrix(x$pstate))
                temp <- x$pstate[igrp==i,, drop=FALSE]
            else if (is.array(x$pstate))
                temp <- matrix(x$pstate[igrp==i,,,drop=FALSE],
                               ncol= nstate)
            else temp <- matrix(x$pstate[igrp==i], ncol=1)

            tt <- x$time[igrp==i]
 
            # Now cut it off at maxtime
            delta <- diff(c(tt[tt<maxtime[i]], maxtime[i]))
            if (length(delta) > nrow(temp)) delta <- delta[1:nrow(temp)]
            if (length(delta) < nrow(temp))
                delta <- c(delta, rep(0, nrow(temp) - length(delta)))
            meantime[i,] <- colSums(delta*temp)

            if (!is.null(x$influence)) {
                # calculate the variance
                if (is.list(x$influence))
                    itemp <- apply(x$influence[[i]], 1,
                                   function(x) colSums(x*delta))
                else itemp <- apply(x$influence, 1,
                                    function(x) colSums(x*delta))
                stdtime[i,] <- sqrt(rowSums(itemp^2))
           }
        }
        outmat <- cbind(outmat, c(meantime)/scale)
        cname <- c("n", "nevent", "rmean")
        if (!is.null(x$influence)) {
            outmat <- cbind(outmat, c(stdtime)/scale)
            cname <- c(cname, "std(rmean)")
        }
        # report back a single time, if there is only one
        if (all(maxtime == maxtime[1])) maxtime <- maxtime[1]
    }
    else cname <- c("n", "nevent")
    dimnames(outmat) <- list(rowname, cname)

    if (rmean=='none') list(matrix=outmat)
    else list(matrix=outmat, end.time=maxtime/scale)
}
"[.survfitms" <- function(x, ..., drop=FALSE) {
    nmatch <- function(i, target) { 
        # This function lets R worry about character, negative, 
        # or logical subscripts
        #  It always returns a set of positive integer indices
        temp <- seq(along.with=target)
        names(temp) <- target
        temp[i]
    }
    if (!is.null(x$influence.pstate) || !is.null(x$influence.cumhaz))
        x <- survfit0(x, x$start.time)  # make influence and pstate align
    ndots <- ...length()      # the simplest, but not avail in R 3.4
    # ndots <- length(list(...))# fails if any are missing, e.g. fit[,2]
    # ndots <- if (missing(drop)) nargs()-1 else nargs()-2  # a workaround
 
    dd <- dim(x)
    dmatch <- match(c("strata", "data", "states"), names(dd), nomatch=0)
    if (is.null(x$states)) stop("survfitms object has no states component")
    if (dmatch[3]==0) stop ("survfitms object has no states dimension")
    dtype <- match(names(dd), c("strata", "data", "states"))

    if (ndots==0) return(x)  # no subscript given
    if (ndots >0 && !missing(..1)) i <- ..1 else i <- NULL
    if (ndots> 1 && !missing(..2)) j <- ..2 else j <- NULL
    if (ndots> 2 && !missing(..3)) k <- ..3 else k <- NULL
    if (is.null(i) & is.null(j) & is.null(k)) return(x) # only one curve
    
    # Make a new object
    newx <- vector("list", length(x))
    names(newx) <- names(x)
    for (kk in c("logse", "version", "conf.int", "conf.type", "type", 
                 "start.time", "call"))
        if (!is.null(x[[kk]])) newx[[kk]] <- x[[kk]]
    newx$transitions <- NULL # may no longer be accurate, and not needed
    class(newx) <- class(x)

    # Like a matrix, let the user use a single subscript if they desire
    if (ndots==1 && length(dd) > 1) {
        # the 'treat it as a vector' case
        if (!is.numeric(i))
            stop("single subscript must be numeric")
        if (any(dmatch==2)) stop("single index subscripts are not supported for a survfit objet with both data and state dimesions")

        # when subscripting a mix, these don't endure
        newx$cumhaz <- newx$std.chaz <- newx$influence.chaz <- NULL
        newx$transitions <- newx$states <- newx$newdata <- NULL
        
        # what strata and columns do I need?
        itemp <- matrix(1:prod(dd), nrow=dd[1])
        jj <- (col(itemp))[i]    # columns
        ii <- (row(itemp))[i]    # this is now the strata id
        
        if (dtype[1]!=1 || dd[1]==1) # no strata or only 1
            irow <- rep(seq(along.with= x$time), length(ii))
        else {
            itemp2 <- split(1:sum(x$strata), rep(1:length(x$strata), x$strata))
            irow <- unlist(itemp2[ii])  # rows of the pstate object
        }
        inum <- x$strata[ii]        # number of rows in each ii
        indx <- cbind(irow, rep(jj,ii))      # matrix index for pstate
        
        # The n.risk, n.event, .. matrices dont have a newdata dimension.
        if (all(dtype!=2) || dd["data"]==1) kk <- jj
        else {  # both data and states
            itemp <- matrix(1:(dd["data"]*dd["states"]), nrow=dd[2])
            kk <- (col(itemp))[jj]    # the state of each selected one
            indx2 <- cbind(irow, rep(k, irow))  
        }
        newx$n <- x$n[ii]
        newx$time <- x$time[irow]
        for (z in c("n.risk", "n.event", "n.censor", "n.enter"))
            if (!is.null(x[[z]])) newx[[z]] <- (x[[z]])[indx2]
        for (z in c("pstate", "std.err", "upper", "lower"))
            if (!is.null(x[[z]])) newx[[z]] <- (x[[z]])[indx]
        
        newx$strata <- x$strata[ii]
        names(newx$strata) <- seq(along.with=ii)
        
        return(newx)
    }
        
    # not a single subscript, i.e., the usual case
    # Backwards compatability: If x$strata=NULL, it is a semantic argument
    #  of whether there is still "1 stratum".  I have used the second
    #  form at times, e.g.  x[1,,2] for an object with only data and state
    #  dimensions.
    # If there are no strata, 1 too many subscripts, and the first is 1,
    #  assume this case and toss the first
    if (ndots == (length(dd)+1)) {
        if (is.null(x$strata) && (is.null(i) || (length(i)==1 && i==1))) {
            i <-j; j <-k; k <- NULL
        } else stop("incorrect number of dimensions")
    } else if (ndots != length(dd)) stop("incorrect number of dimensions")

    # create irow, which selects for the time dimension of x
    if (dtype[1]!=1 || is.null(i)) {
        irow <- seq(along.with= x$time)
    }
    else {
        i <- nmatch(i, names(x$strata))
        itemp <- split(1:sum(x$strata), rep(1:length(x$strata), x$strata))
        irow <- unlist(itemp[i])  # rows of the pstate object
    }

    # Select the n, strata, and time components of the output.  Make j,k
    #  point to the subscripts other than strata (makes later code a touch
    #  simpler.)
    newx$time <- x$time[irow]    
    if (dtype[1] !=1) {  # there are no strata
        newx$n <- x$n
        k <- j; j <- i;
        dd <- c(0, dd)
        dtype <- c(1, dtype)
    }
    else { # there are strata
        if (is.null(i)) i <-seq(along.with=x$strata)
        if ((drop && length(i)>1) || !drop) newx$strata <- x$strata[i]
        newx$n <- x$n[i]
    }

    # The n.censor and n.enter values do not repeat with multiple X values
    for (z in c("n.censor", "n.enter"))
        if (!is.null(x[[z]])) newx[[z]] <- (x[[z]])[irow, drop=FALSE]
    
    # two cases: with newx or without newx  (pstate is always present)
    nstate <- length(x$states)
    if (dtype[2] !=2) {  # j indexes the states, there is no data dimension
        if (is.null(j)) j <- seq.int(nstate)
        else j <- nmatch(j, x$states)

        # keep these as start points for plotting, even though they won't make
        #  true sense if states are subset, since rows won't sum to 1
        if (!is.null(x$p0)) {
            if (is.matrix(x$p0)) newx$p0 <- x$p0[i,j, drop=FALSE] 
            else newx$p0 <- x$p0[j]
        }
        if (!is.null(x$sp0)) {
            if (is.matrix(x$sp0)) newx$sp0 <- x$sp0[i,j, drop=FALSE] 
            else newx$sp0 <- x$sp0[j]
        }   
        
        # in the rare case of a single strata with 1 obs, don't drop dims
        if (length(irow)==1 && length(j) > 1) drop2 <- FALSE 
        else drop2 <- drop

        for (z in c("n.risk", "n.event"))
            if (!is.null(x[[z]])) newx[[z]] <- (x[[z]])[irow,j, drop=drop2]
        for (z in c("pstate", "std.err", "upper", "lower"))
            if (!is.null(x[[z]])) newx[[z]] <- (x[[z]])[irow,j, drop=drop2]
        if (!is.null(x$influence.pstate)) {
            if (is.list(x$influence.pstate)) {
                if (length(i)==1) newx$influence.pstate <- x$influence.pstate[[i]]
                else newx$influence.pstate <- lapply(x$influence.pstate[i],
                                     function(x) x[,,j, drop= drop])
                }
            else newx$influence.pstate <- x$influence.pstate[,,j, drop=drop]
        }

        if (length(j)== nstate && all(j == seq.int(nstate))) {
            # user kept all the states, in original order
            newx$states <- x$states
            for (z in c("cumhaz", "std.chaz"))
                 if (!is.null(x[[z]])) newx[[z]] <- (x[[z]])[irow,, drop=drop2]
            if (!is.null(x$influence.chaz)) {
                if (is.list(x$influence.chaz)) {
                    newx$influence.chaz <- x$influence.chaz[i]
                    if (length(i)==1 && drop) 
                        newx$influence.chaz <- x$influence.chaz[[i]]
                }
                else newx$influence.chaz <- x$influence.chaz
            }
        }
        else {
            # Some states were dropped, leaving no consistent way to 
            #  subscript cumhaz, or not one I have yet seen clearly
            # So remove it from the object
            newx$cumhaz <- newx$std.chaz <- newx$influence.chaz <- NULL
            if (length(j)==1 & drop) {
                newx$states <- NULL
                temp <- class(newx)
                class(newx) <- temp[temp!="survfitms"]
            }
            else newx$states <- x$states[j]
        }
    }
    else {  # j points at newdata, k points at states
        if (is.null(j)) j <- seq.int(dd[2])
        else j <- nmatch(j, seq.int(dd[2]))

        if (is.null(k)) k <- seq.int(nstate)
        else k <- nmatch(k, x$states)

        # keep these as start points for plotting, even though they won't make
        #  true sense is states are subset, since rows won't sum to 1
        # (all data= sets have the same p0)
        if (!is.null(x$p0)) {
            if (is.matrix(x$p0)) newx$p0 <- x$p0[i,k] else newx$p0 <- x$p0[k]
        }
        if (!is.null(x$sp0)) {
            if (is.matrix(x$sp0)) newx$sp0 <- x$p0[i,k] else newx$sp0 <- x$sp0[k]
        }   

         if (length(irow)==1) {
            if (length(j) > 1) drop2 <- FALSE else drop2<- drop
            if (length(k) > 1) drop3 <- FALSE else drop3 <- drop
        } 
        else drop2 <- drop3 <- drop

        for (z in c("n.risk", "n.event"))
            if (!is.null(x[[z]])) newx[[z]] <- (x[[z]])[irow, k, drop=drop3]
        for (z in c("pstate", "std.err", "upper", "lower"))
            if (!is.null(x[[z]])) newx[[z]] <- (x[[z]])[irow,j,k, drop=drop2]
  
        if (!is.null(x$influence.pstate)) {
            if (is.list(x$influence.pstate)) {
                if (length(i)==1) 
                    newx$influence.pstate <- (x$influence.pstate[[i]])[,,j,k, drop=drop]
                else newx$influence.pstate <- lapply(x$influence.pstate[i],
                                     function(x) x[,,j,k, drop= drop])
                }
            else newx$influence.pstate <- x$influence.pstate[,,j,k, drop=drop]
        }

        if (length(k)== nstate && all(k == seq.int(nstate))) {
            # user kept all the states
            newx$states <- x$states
            for (z in c("cumhaz", "std.chaz"))
                 if (!is.null(x[[z]])) 
                     newx[[z]] <- (x[[z]])[irow,j,, drop=drop2]
            if (!is.null(x$influence.chaz)) {
                if (is.list(x$influence.chaz)) {
                    newx$influence.chaz <- (x$influence.chaz[i])[,j,]
                    if (length(i)==1 && drop) 
                        newx$influence.chaz <- x$influence.chaz[[i]]
                }
                else newx$influence.chaz <- x$influence.chaz[,j,]
            }
        }
        else {
            # never drop the states component.  Otherwise downstream code
            #  will start looking for x$surv instead of x$pstate
            newx$states <- x$states[k]
            newx$cumhaz <- newx$std.chaz <- newx$influence.chaz <- NULL
            x$transitions <- NULL
         }

        if (length(j)==1 && drop) newx$newdata <- NULL
        else newx$newdata <- x$newdata[j,,drop=FALSE]  #newdata is a data frame
 
    }
    newx
}
