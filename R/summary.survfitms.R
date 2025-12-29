# Summary function for survfitms and survfit.coxphms objects
# Very similar to summary.survfit: substitute pstate for surv, and more objects
#  are matrices
summary.survfitms <- function(object, times, censored=FALSE, 
                            scale=1, extend=FALSE, 
                            rmean=getOption('survfit.rmean'),
                            data.frame= FALSE, dosum=TRUE,
                            ...) {
    if (!inherits(object, 'survfit'))
            stop("summary.survfit can only be used for survfit",
                 " and survfit.coxph objects")
    fit <- object  # I get tired of typing "object"
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

    fit0 <- survfit0(object)  #add time 0
    if (!data.frame) {
        # adding time 0 makes the mean and median easier
        temp <- survmean2(fit0, scale=scale, rmean)  
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
    #  and is not summed, to pick of the right rows of times
    # the pmax(1, indx) is for when a user gives a time before the first event
    #  in which case findInterval has given a 0
    ssub<- function(x, indx) {  #select an object and index
        if (is.logical(indx)) indx <- which(indx)
        if (!is.null(x) && length(indx)>0) {
            if (is.matrix(x)) x[pmax(1,indx),,drop=FALSE]
            else x[pmax(1, indx)]
        }
        else NULL
    }
    
    # By replacing elements of fit, we retain the ones that don't need
    #  to be subscripted (n, logse, etc)
    if (missing(times)) {
        if (!censored) {  # do not retain time points with no events
            index <- (rowSums(fit$n.event) >0)
            for (i in c("time","n.risk", "n.event", "pstate", "std.err", 
                                "upper", "lower", "cumhaz", "std.chaz")) {
                if (!is.null(fit[[i]])) {  # not all components in all objects
                    fit[[i]] <- ssub(fit[[i]], index)
                }
            }

            # The n.enter, n.censor, and n.transition values are accumualated
            #  all are matrices
            if (is.null(fit$strata)) {
                for (i in c("n.enter", "n.censor", "n.transition"))
                    if (!is.null(fit[[i]]))
                        fit[[i]] <- delta(fit[[i]], index)
            }
            else {
                sindx <- rep(1:nstrat, fit$strata)
                for (i in c("n.enter", "n.censor", "n.transition")) {
                    if (!is.null(fit[[i]]))
                        fit[[i]] <- unlist(lapply(1:nstrat, function(j) 
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
    else { # the user specified times
        if (length(times) ==0) stop("no values in times vector")
        if (inherits(times, "Date")) times <- as.numeric(times) # allow Dates
        if (!is.numeric(times)) stop("times must be a numeric vector")
        if (!all(is.finite(times))) stop("times contains missing or infinite values")
        times <- unique(sort(times))
        fit <- fit0  # makes summing easier 

        # findrow is called once per stratum
        # We have to deal with user specified times that are before the first
        #  time value in the curve or after the last, which is easier done one
        #  curve at at time.  Ditto for summing number of events & censors.
        # For the summed terms (n.event, n.censor, n.enter) it turns out to
        #  be easier to extract them within findrow, ditto for the updated
        #  time vector (a user might ask for times that are past the end of
        #  one stratum, but not past the end of another).  
        # For other components it is easier to return the index vectors 
        #  do the subset after findrow, all at once.
        # The calling routine makes use of the fact that [.survfit knows how
        #  to extract a single stratum, giving a simpler curve.
        findrow <- function(fit, times, extend, first) {
            if (!extend) {
                maxtime <- max(fit$time)
                times <- times[times <= maxtime]
            }
            ntime <- length(fit$time)
            if (ntime ==0) {
                return(list(time=NULL, index1 =NULL, index2=NULL, 
                            n.censor=NULL,  n.event=NULL, n.enter=NULL))
            }
                            
            index1 <- findInterval(times, fit$time) 
            index2 <- 1 + findInterval(times, fit$time, left.open=TRUE)
            index2 <- ifelse(index2 > length(fit$time), 0, index2 + first)  

            flist <- list(time =times, index1= index1+first, index2= index2)
            # Now the ones that are summed between reporting times
            for (i in c("n.event", "n.censor", "n.enter", "n.transition")) {
                if (dosum) flist[[i]] <- delta(fit[[i]], index1)
                else flist[[i]] <- ssub(fit[[i]], index1)
            }
            # within this function index1 = index to the rows of this stratum
            # returned value = index to the rows of the full survfit object
            flist
        }

        if (nstrat ==1) {
            ltemp <- list(findrow(fit, times, extend, 0))
        }
        else {
            # first = first obs of each stratum
            first <- cumsum(c(0, fit$strata[-length(fit$strata)]))
            ltemp <- vector("list", nstrat)
            if (length(dim(fit)) > 2) { # there is a data dimension
                for (i in 1:nstrat) 
                    ltemp[[i]] <- findrow(fit[i,,], times, extend, first[i])
            } else { 
                for (i in 1:nstrat) 
                    ltemp[[i]] <- findrow(fit[i,], times, extend, first[i])
            }
        }
 
        # Replace bits of the fit object one by one
        # get() saves me some typing
        # For strata, if the times arg ended up with 0 rows in a stratum,
        #  we retain that strata with a count of 0. This is resolved further
        #  below when strata is expanded.
        get <- function(x,y) lapply(x, function(x) x[[y]]) 
        if (!is.null(fit$strata))
            fit$strata[] <- sapply(get(ltemp, "time"), length)
        fit$time <- unlist(get(ltemp, "time"))
        
        # Elements n, n.id, p0, logse, conf.type, conf.ing, states, type, t0,
        #  call are left as is.
        # The summed elements need to be unstacked
        for (i in c("n.event", "n.censor", "n.enter", "n.transition")) {
            if (!is.null(fit[[i]]))
                fit[[i]] <- do.call(rbind, get(ltemp, i))
        }

        # For an intermediate time, e.g., curve has time points of 10 and 20,
        #  user asked for 15, the n.risk component will be for the later
        #  time.  Every curve ends with censors or deaths, so n.risk beyond
        #  the endpoint is 0
        fit$n.risk <- rbind(0, fit$n.risk)[1L + unlist(get(ltemp, "index2")),]

        # For the curves, an intermediate time maps to the earlier
        # time point, e.g., the pstate at time 15 is the curve value at time 10
        index1 <- unlist(get(ltemp, "index1"))
        for (i in c("pstate", "upper", "lower", "std.err", "cumhaz",
                    "std.chaz")) {
            if (!is.null(fit[[i]])) {
                # matrix= survfitms, array= survfit.coxphms
                if (is.matrix(fit[[i]]))
                    fit[[i]] <-(fit[[i]])[index1,, drop=FALSE]
                else  fit[[i]] <-(fit[[i]])[index1,,, drop=FALSE]
            }
        }
    }

    # finish off the output structure
    if (scale != 1) {
        # fix scale in the output
        fit$time <- fit$time/scale
    }
    # Strata is expanded to an element per row
    if (!is.null(fit$strata)) 
        fit$strata <- factor(rep(1:nstrat, fit$strata), 1:nstrat,
                             labels= names(fit$strata))
    if (data.frame) {
        # dim(fit) will be (strata, data, state), the data frame should be
        #  in the standard 'first subscript varies fastest' order used by R.
        #  Think of (time,strata) as a single subscript.
        # The cumulative hazard, std.chaz, and n.transitions component have
        #  a different dimension (number of transitions vs number of states)
        #  so don't get included, and likewise informational things like call,
        #  t0, n, n.id, p0, ...
        # The array bits of the fit (pstate, std.err) can be strung out with c()
        # The matrix ones (n.risk, n.event, n.censor) have dimensions of
        #  (time-strata, state) so need to be expanded by making nrow(newdata)
        #  copies of each column.
        # The newdata expansion is the most complex, e
        if (is.null(fit$newdata)) nd <-1
        else  nd <- nrow(fit$newdata) # data dimension, might be 1
        nt <- nrow(fit$n.risk) # the time-strata dimension
        ns <- length(fit$states)
        
        j <- rep(1:nt, nd)
        new <- data.frame(time= rep(fit$time, nd*ns),
                          n.risk  = c(fit$n.risk[j,]),
                          n.event = c(fit$n.event[j,]),
                          n.censor= c(fit$n.censor[j,]),
                          pstate = c(fit$pstate))
        if (!is.null(fit$std.err)) new <- cbind(new, std.err= c(fit$std.err),
                          lower = c(fit$lower), upper = c(fit$upper))
        if (!is.null(fit$strata)) new$stratra <- rep(fit$strata, nd*ns)
        new$state <- rep(fit$states, each= nd*nt)
        if (!is.null(fit$newdata)) { #coxph curves
            k <- rep(rep(1:nd, each=nt), ns)
            new <- cbind(new, fit$newdata[k,])
        }
        new  # result is of class data.frame
    } else {
        fit$table <- table
        if (length(rmean.endtime)>0  && !any(is.na(rmean.endtime[1]))) 
            fit$rmean.endtime <- rmean.endtime
        class(fit) <- "summary.survfitms"
        fit
    }
}
