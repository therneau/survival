#
# Get pseudo values for a survival curve, based on the IJ
#
pseudo <- function(fit, times, type, collapse=TRUE, data.frame=FALSE, ...){
    if (!inherits(fit, "survfit"))
        stop("fit argument must be a survfit object")
    # Add the model.frame here, rather than one step down in the residuals fcn
    # It seems to help some "not found" type errors.
    if (is.null(fit$model)) # add it before calling residuals
        fit$model <- model.frame(fit)

    if (missing(type)) type <- "pstate"
    else {
        # allow various aliases for the type
        temp <- c("pstate", "cumhaz", "sojourn", "survival",
                              "chaz", "rmst", "rmts", "auc")
        type <- match.arg(casefold(type), temp)
        itemp <-  c(1,2,3,1,2,3,3,3)[match(type, temp)]
        type <- c("pstate", "cumhaz", "auc")[itemp]
    }
    ntime <- length(times)

    # all the real work is done by the residuals function
    rfit <- residuals.survfit(fit, times=times, type=type, array=FALSE,
                              collapse= collapse, data.frame=data.frame, 
                              extra= TRUE, ...)
    # The residuals are modifed by multiplying by the number of subjects
    # and then centering at that survival estimate.
    # Set up the indexing so that the right survival estimate is used for
    # each residual
 
    ddf <- dim(fit)
    if (!data.frame) ddr <- dim(rfit$resid)
    if (is.null(fit$strata)) {
        ncurve <- 1L  
        icurve <- 1L
        
    } else {
        ncurve <- length(fit$strata)  # this should match ddf[1]
        icurve <- rfit$curve
    }
    
    if (!is.null(fit$states)) {
        nstate <- length(fit$states)
        nhaz <- ncol(fit$cumhaz)
        if (data.frame) {
            if (type == "cumhaz") 
                ihaz <- match(rfit$transition, colnames(fit$cumhaz))
            else istate <- match(rfit$state, fit$states)
        } else {
            # the residuals are a matrix with subject, state or transition,
            #  and times as dimensions
            if (type== "cumhaz") ihaz <- rep(1:nhaz  , each= ddr[1])
            else  istate <- rep(1:nstate, each= ddr[1])
        }
    } else {
        nstate <- nhaz <- 1L
        istate <- ihaz <- 1L
    }

    if (ntime ==1)  itime <- 1L
    else {
        if (data.frame) itime <- match(rfit$time, times)   
        else itime <- rep(1:ntime, each= prod(ddr[-length(ddr)]))
    }

    # The residuals are inflated by nn
    if (is.null(fit$n.id)) nn <- fit$n  # each row a separate subject
    else nn <- fit$n.id
   
    # warn about pseudo values after the end of a curve
    if (ncurve==1) etime <- max(fit$time)
    else etime <- fit$time[cumsum(fit$strata)]
    if (any(etime < max(times))) 
        warning("requested time points are beyond the end of one or more curves")

    # Now do the work
    if (type == "pstate") { 
        temp <- summary(fit, times=times, extend=TRUE)
        if (nstate >1) { # multi-state
            # pstate will have times, curve, state,  as the order
            # resid will have subject, state, times
            yhat <-  array(temp$pstate, dim=c(ntime, ncurve, nstate))
            ptemp <- yhat[cbind(itime, icurve, istate)] + nn[icurve]* rfit$resid
        }
        else {
            # yhat will be a matrix with one row per time, one col per curve
            yhat <- matrix(temp$surv, ncol = ncurve)
            ptemp <- yhat[cbind(itime, icurve)] + nn[icurve]*rfit$resid
        }
    } else if (type =="auc") { 
        if (nstate >1) { # multi-state
            yhat <- array(0, dim=c(ntime, ncurve, nstate))
            # summary.survfit gives only one AUC time per call
            for (i in 1:ntime) {
                temp <- summary(fit, rmean=times[i])$table
                yhat[i,,] <- temp[,"rmean"]
            }
            ptemp <- yhat[cbind(itime, icurve, istate)] + nn[icurve]* rfit$resid
        }
        else {  
            # yhat will be a matrix with one row per time, one col per curve
            yhat <- matrix(0, ntime, ncurve)
            for (i in 1:ntime) {
                temp <- summary(fit, rmean=times[i])$table 
                if (ncurve==1) yhat[i,] <- temp["rmean"]
                else yhat[i,] <- temp[,"rmean"]
            }
            ptemp <- yhat[cbind(itime, icurve)] + nn[icurve]*rfit$resid
        }
    } else if (type== "cumhaz") {
        temp <- summary(fit, times=times, extend=TRUE)
        if (nstate > 1) {
            yhat <- array(temp$cumhaz, dim=c(ntime, ncurve, nhaz))
            ptemp <- yhat[cbind(itime, icurve, ihaz)] + nn[icurve]* rfit$resid
        }
         else {
            yhat <- matrix(temp$cumhaz, ncol=ncurve)
            ptemp <- yhat[cbind(itime, icurve)] + nn[icurve]* rfit$resid
        }
    } else stop("unknown type")  # this should never happen

    if (data.frame) { 
        rfit$pseudo <- ptemp
        rfit
    } else drop(ptemp)
}       

