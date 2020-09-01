# Function to compute Brier scores

brier <- function(fit, times, newdata, censformula) {
    if (!inherits(fit, "coxph"))
        stop("fit argument must be a coxph model")
    if (inherits(fit, "coxphms"))
        stop("not yet applicable to multi-state models")

    yform <- formula(dfit1)  # the formula in the coxph object
    ny <- ncol(fit$y)

    if (!missing(censformula)) {
        if (length(censformula) <3) 
            stop("censoring formula must have a response")
        if (missing(newdata))
            stop("a censoring formula requires newdata")
    }

    if (missing(newdata)) {
        # Baseline predicted probabilities
        s0 <- survfit(fit$y ~1, se.fit=FALSE)
        s1 <- survfit(fit, newdata= model.frame(fit), se.fit=FALSE)  

        ytemp <- fit$y
        ytemp[,ny] <- 1-ytemp[,ny]
        c0 <- survfit0(survfit(yemp ~1))
        dstat <- fit$y[,ny]
    }
    else {  # apply to new data
        tform <- yform
        tform[[3]] <- 1
        s0 <- survfit(yform, data=newdata, se=FALSE)
        s1 <- survfit(fit, newdata=newdata, se.fit=FALSE)
        if (missing(censformula))
            stop("not yet written")
        else c0 <- survfit0(survfit(censformula, data=newdata))
        dstat <- eval(tform[[1]], newdata)
    }

    if (missing(times)) times <- s0$time[s0$n.event >0]
    p0 <- 1- summary(s0, times)$surv
    p1 <- 1- summary(s1, times)$surv

    n <- length(dstat)
    ntime <- length(times)
    # Censoring distribution
    cx <- c0$time
    cy <- c0$surv

    b0 <- b1 <- matrix(0, nrow=n, ncol=ntime)
    wt <- b0
    for (i in 1:ntime) {
        b0[,i] <- ifelse(dtime > times[i], p0[i]^2, (dstat- p0[i])^2)
        b1[,i] <- ifelse(dtime > times[i], p1[i,]^2, (dstat- p1[i,])^2)
        j <- findInterval(pmin(dtime, times[i]), cx, left.open=TRUE)
        wt[,i] <- ifelse(dtime< times[i] & dstat==0, 0, 1/cy[j])
    }       

    brier <- cbind(colMeans(b0*wt),
                   colMeans(b1*wt))
    dimnames(brier) <- list(times, c("NULL", "Model"))
    return(list(brier=brier, n= colSums(wt>0), p0 = p0, times=times,
                b0= b0, b1=b1, wt=wt))
}
