# Summary function for survfit and survfit.coxph objects
summary.survfit <- function(object, times, censored=FALSE, 
                            scale=1, extend=FALSE, 
                            rmean=getOption('survfit.rmean'),
                            data.frame= FALSE, 
                            ...) {
    fit <- object  # I get tired of typing "object"
    if (!inherits(fit, 'survfit'))
            stop("summary.survfit can only be used for survfit",
                 " and survfit.coxph objects")
    if (is.null(fit$logse)) fit$logse <- TRUE   #older style objects lack this

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

    fit0 <- survfit0(fit)
    if (!data.frame) {
        # adding time 0 makes the mean and median easier
        temp <- survmean(fit0, scale=scale, rmean)  
        table <- temp$matrix  #for inclusion in the output list
        rmean.endtime <- temp$end.time
    }

    if (!is.null(fit$strata)) {
        nstrat <-  length(fit$strata)
    } else nstrat <- 1

    # If times is present, then n.event, n.censor, and n.enter are summed
    #  between those time points.  Utility function to do that
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

    # called for each component of the curve that has a time dimension
    #  and is not summed
    ssub<- function(x, indx) {  #select an object and index
        if (is.logical(indx)) indx <- which(indx)
        if (!is.null(x) && length(indx)>0) {
            if (is.matrix(x)) x[pmax(1,indx),,drop=FALSE]
            else if (is.array(x))  x[pmax(1,indx),,,drop=FALSE]
            else x[pmax(1, indx)]
        }
        else NULL
    }
    
    # By replacing components of fit, summary.surfit inherits several bits
    if (missing(times)) {
        if (!censored) {  # do not retain time points with no events
            index <- fit$n.event >0
            for (i in c("time","n.risk", "n.event", "surv", "std.err", 
                                "upper", "lower", "cumhaz", "std.chaz")) {
                if (!is.null(fit[[i]])) {  # not all components in all objects
                    fit[[i]] <- ssub(fit[[i]], index)
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
                #  events at all, only censored and hence 0 lines of output
                # the [] retains the original names
                fit$strata[] <- as.vector(table(factor(sindx[index], 1:nstrat))) 
            }
        }
        #if missing(times) and censored=TRUE, the fit object is ok as it is
    }
    else {
        if (length(times) ==0) stop("no values in times vector")
        if (inherits(times, "Date")) times <- as.numeric(times) # allow Dates
        if (!is.numeric(times)) stop("times must be a numeric vector")
        if (!all(is.finite(times))) stop("times contains missing or infinite values")  
        times <- unique(sort(times))
        fit <- fit0  # findrow() needs the starting time

        # findrow is called once per stratum
        #   times will be the user specified times
        #   returned is a subset of the rows for the stratum
        # We have to deal with user specified times that are before the first
        #  time value in the curve or after the last, which is easier done one
        #  curve at at time
        findrow <- function(fit, times, extend) {
            if (!extend) {
                maxtime <- max(fit$time)
                times <- times[times <= maxtime]
            }
            ntime <- length(fit$time)
            if (ntime ==0) { 
                if (data.frame) return(list(time = times))
                else stop("no points selected for one or more curves,", 
                     " data error (?) or consider using the extend argument")
            }
                            
            index1 <- findInterval(times, fit$time) 
            index2 <- 1 + findInterval(times, fit$time, left.open=TRUE)
                
            fit$time <- times
            for (i in c("surv", "upper", "lower", "std.err", "cumhaz",
                        "std.chaz")) {
                if (!is.null(fit[[i]])) fit[[i]] <- ssub(fit[[i]], index1)
            }
            
            # Every observation in the data has to end with a censor or event.
            #  So by definition the number at risk after the last observed time
            #  value must be 0.
            fit$n.risk <- c(fit$n.risk, 0)[index2]

            for (i in c("n.event", "n.censor", "n.enter"))
                fit[[i]] <- delta(fit[[i]], index1)
            fit
        }

        if (nstrat ==1) fit <- findrow(fit, times, extend)
        else {
            ltemp <- vector("list", nstrat)
            if (length(dim(fit)) > 1) {
                for (i in 1:nstrat) 
                    ltemp[[i]] <- findrow(fit[i,], times, extend)
            } else { 
                for (i in 1:nstrat) ltemp[[i]] <- findrow(fit[i], times, extend)
            }
         
            # now stack them: time= c(time for curve 1, time for curve 2, etc)
            #  and so on for all components
            unlistsurv <- function(x, name) {
                temp <- lapply(x, function(x) x[[name]])
                if (is.vector(temp[[1]])) unlist(temp)
                else if (is.matrix(temp[[1]])) do.call("rbind", temp)
            }

            # unlist all the components built by a set of calls to findrow
            #  and remake the strata
            keep <- c("time", "surv", "upper", "lower", "std.err",
                      "cumhaz", "n.risk", "n.event", "n.censor", "n.enter",
                      "std.chaz")
            for (i in keep) 
                if (!is.null(fit[[i]])) fit[[i]] <- unlistsurv(ltemp, i)
            fit$strata[] <- sapply(ltemp, function(x) length(x$time))
        }
    }

    # finish off the output structure
    # A survfit object may contain std(log S) or std(S), summary always std(S)
    if (!is.null(fit$std.err) && fit$logse) 
        fit$std.err <- fit$std.err * fit$surv   
    if (scale != 1) {
        # fix scale in the output
        fit$time <- fit$time/scale
    }

    if (data.frame) {
        fit <- unclass(fit)  # toss the survfit class
        indx <- match(c("time", "n.risk", "n.event", "n.censor", 
                        "surv", "cumhaz", "std.err", "std.chaz",
                        "lower", "upper"), names(fit), nomatch=0)
        if (!is.null(fit$strata))
            newstrat <- factor(rep(1:nstrat, fit$strata), 1:nstrat,
                               labels= names(fit$strata))
        if (is.matrix(fit$surv)) { # survfit.coxph object
            nc <- ncol(fit$surv)
            ndata <- lapply(fit[indx], function(x) {
                                 if (length(x)==0) NULL
                                 else if (is.matrix(x)) c(x)
                                 else rep(x, nc)})
            ndata <- data.frame(ndata)
            if (!is.null(fit$strata)) 
                ndata$strata <- rep(newstrat, nc)
            ndata$data <- rep(1:nc, each= length(fit$time))
        } else {
            ndata <- data.frame(fit[indx])
            if (!is.null(fit$strata)) ndata$strata <- newstrat
        }
        ndata
    } else {
        fit$table <- table
        if (length(rmean.endtime)>0  && !any(is.na(rmean.endtime[1]))) 
            fit$rmean.endtime <- rmean.endtime
        # Expand the strata. It has used 1,2,3 for a long while
        if (!is.null(fit$strata)) 
            fit$strata <- factor(rep(1:nstrat, fit$strata), 1:nstrat,
                                 labels= names(fit$strata))
        class(fit) <- "summary.survfit"
        fit
    }
}
