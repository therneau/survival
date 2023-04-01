brier <- function(fit, times, newdata, ties=TRUE, detail =FALSE, timefix=TRUE,
                  efron = FALSE) {
    Call <- match.call()
    if (!inherits(fit, "coxph")) stop("fit must be a coxph object")

    if (missing(newdata)) mf <- stats::model.frame(fit)
    else mf <-stats::model.frame(fit, data=newdata)
    Y <- mf[[1]] # the survival object
    if (!is.Surv(Y)) stop("response must be a Surv object")
    type <- attr(Y, "type")
    if (!(type %in% c("right", "mright", "counting", "mcounting")))
        stop("response must be right censored")
    n <- nrow(Y)
    ny <- ncol(Y)  # 3 = time1, time2 data

    if (!is.logical(timefix) || length(timefix) > 1)
        stop("invalid value for timefix option")
    if (timefix) Y <- aeqSurv(Y) 

    casewt <- model.weights(mf)
    if (is.null(casewt)) casewt <- rep(1, n)
    else {
        if (!is.numeric(casewt)) stop("weights must be numeric")
        if (any(!is.finite(casewt))) stop("weights must be finite") 
        if (any(casewt <0)) stop("weights must be non-negative")
        casewt <- as.numeric(casewt)  # transform integer to numeric
    }

    id <- model.extract(mf, "(id)")
    if (ny==3 && (is.null(id))) stop("id is required for start-stop data")
    if (!is.null(id)) {
        if (is.null(attr(Y, 'states'))) {
            ytemp <- Y
            attr(ytemp, 'states') <- 'event'  # survcheck2 wants a states attr
            check <- survcheck2(ytemp, id)
        }
        else check <- survcheck2(Y, id)

        if (any(check$flag > 0)) 
                stop("one or more flags are >0 in survcheck")
        n.startstate <- sum(check$transitions[,1] >1)
        if (ny ==2) samestart=TRUE   # everyone starts at the same time
        else {
           etemp  <- tapply(Y[,1], id, min)
           samestart <- all(etemp== etemp[1])
        }    
    } else check <- NULL
    if (length(id)==0) id <- seq.int(n)

    # Finally, it's time to do the work.
    if (is.null(check) || (n.startstate==1 & samestart)) simple <- TRUE
    else simple <- FALSE
    if (!simple) stop("delayed entry is not yet implemented")
    
    # For baseline predicted probabilities, use the same approximations as the
    #  Cox model, if allowed by the 'efron' option.
    if (efron && fit$method == "efron") 
        s0 <- survfit(Y ~1, weights= casewt, se.fit=FALSE, ctype=2, stype=2)
    else s0 <- survfit(Y ~1, weights= casewt,se.fit=FALSE, stype=1)
    if (missing(times)) times <- s0$time[s0$n.event >0]
    p0 <- 1- summary(s0, times, extend=TRUE)$surv

    # model predictions
    if (missing(newdata))
        s1 <- survfit(fit, newdata= fit$call$data, se.fit=FALSE)# FALSE is faster
    else s1 <- survfit(fit, newdata= newdata, se.fit=FALSE)
    p1 <- 1- summary(s1, times, extend=TRUE)$surv

    dtime <- Y[,ny-1]  # time and status in the data
    dstat <- Y[,ny]
    n <- length(dtime)
    ntime <- length(times)

    # Censoring distribution
    if (ties) {
        # move the censorings just a bit to the right
        mindiff <- min(diff(sort(unique(dtime))))
        dtime <- dtime + ifelse(dstat==0, mindiff/2, 0)
    }       
                        
    c0 <- survfit0(survfit(Surv(dtime, 1-dstat) ~ 1, weights=casewt))
    b0 <- b1 <- matrix(0, nrow=n, ncol=ntime)
    casewt <- casewt/sum(casewt)
    brier <- matrix(0, ntime, 2)
    eff.n <- double(ntime)  # the effective n
    for (i in 1:ntime) {
        indx <- findInterval(pmin(dtime, times[i]), c0$time) 
        wt <- ifelse(dtime < times[i] & dstat==0, 0, casewt/c0$surv[indx])
        eff.n[i] <- 1/sum(wt^2)  # the sum of the weights is always 1
                     
        b0[,i] <- ifelse(dtime > times[i], p0[i]^2, (dstat- p0[i])^2)
        b1[,i] <- ifelse(dtime > times[i], p1[i,]^2, (dstat- p1[i,])^2)
        brier[i,] <- c(sum(wt*b0[,i]), sum(wt*b1[,i]))/ sum(wt)
    }       
    ret <- list(rsquared= 1- brier[,2]/brier[,1], brier=brier[,2], times=times)
    if (detail) {
        ret$p0 <- p0
        ret$phat <- p1
        ret$eff.n <- eff.n
    }
    ret
}
