brier <- function(fit, times, newdata, ties=TRUE) {
    # Baseline predicted probabilities, use the same approximations as the
    #  Cox model
    if (fit$method == "efron") 
        s0 <- survfit(fit$y ~1, se.fit=FALSE, ctype=2, stype=2)
    else s0 <- survfit(fit$y ~1, se.fit=FALSE, stype=1)
    if (missing(times)) times <- s0$time[s0$n.event >0]
    p0 <- 1- summary(s0, times, extend=TRUE)$surv

    # model predictions
    s1 <- survfit(fit, newdata= fit$call$data, se.fit=FALSE)# FALSE is faster
    p1 <- 1- summary(s1, times, extend=TRUE)$surv

    ny <- ncol(fit$y)
    dtime <- fit$y[,ny-1]  # time and status in the data
    dstat <- fit$y[,ny]
    n <- length(dtime)
    ntime <- length(times)

    # Censoring distribution
    if (ties) {
        # move the censorings just a bit to the right
        mindiff <- min(diff(sort(unique(dtime))))
        dtime <- dtime + ifelse(dstat==0, mindiff/2, 0)
    }       
                        
    c0 <- survfit0(survfit(Surv(dtime, 1-dstat) ~ 1))
    b0 <- b1 <- matrix(0, nrow=n, ncol=ntime)
    wt <- rep(1/n, n)  # everyone starts out with a weight of 1/n
    brier <- matrix(0, ntime, 2)
    for (i in 1:ntime) {
        indx <- findInterval(pmin(dtime, times[i]), c0$time) 
        wt <- ifelse(dtime < times[i] & dstat==0, 0, 1/c0$surv[indx])
                     
        b0[,i] <- ifelse(dtime > times[i], p0[i]^2, (dstat- p0[i])^2)
        b1[,i] <- ifelse(dtime > times[i], p1[i,]^2, (dstat- p1[i,])^2)
        brier[i,] <- c(sum(wt*b0[,i]), sum(wt*b1[,i]))/ sum(wt)
    }       

    dimnames(brier) <- list(times, c("NULL", "Model"))
    list(brier=brier, b0= b0, b1=b1, time=times)
}
