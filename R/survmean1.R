survmean1 <- function (x, scale = 1, rmean, s = 0.5) 
{
    if (!is.null(x$start.time)) 
        start.time <- x$start.time
    else start.time <- min(0, x$time)
    pfun <- function(nused, time, surv, n.risk, n.event, lower, 
        upper, start.time, end.time) {
        minmin <- function(y, x) {
            keep <- (!is.na(y) & y <= s)
            if (!any(keep)) 
                NA
            else {
                x <- x[keep]
                y <- y[keep]
                if (y[1] == s && any(y < s)) 
                  (x[1] + x[min(which(y < s))])/2
                else x[1]
            }
        }
        if (!is.na(end.time)) {
            hh <- ifelse((n.risk - n.event) == 0, 0, n.event/(n.risk * 
                (n.risk - n.event)))
            keep <- which(time <= end.time)
            temptime <- c(time[keep], end.time)
            tempsurv <- c(surv[keep], surv[max(keep)])
            hh <- c(hh[keep], 0)
            n <- length(temptime)
            delta <- diff(c(start.time, temptime))
            rectangles <- delta * c(1, tempsurv[-n])
            varmean <- sum(cumsum(rev(rectangles[-1]))^2 * rev(hh)[-1])
            mean <- sum(rectangles) + start.time
        }
        else {
            mean <- 0
            varmean <- 0
        }
        med <- minmin(surv, time)
        if (!is.null(upper)) {
            upper <- minmin(upper, time)
            lower <- minmin(lower, time)
            c(nused, max(n.risk), n.risk[1], sum(n.event), sum(mean), 
                sqrt(varmean), med, lower, upper)
        }
        else c(nused, max(n.risk), n.risk[1], sum(n.event), sum(mean), 
            sqrt(varmean), med, 0, 0)
    }
    stime <- x$time/scale
    if (is.numeric(rmean)) 
        rmean <- rmean/scale
    surv <- x$surv
    plab <- c("records", "n.max", "n.start", "events", "*rmean", 
        "*se(rmean)", paste('P',(1-s)*100,sep=''), paste(x$conf.int,
c("LCL", "UCL"), 
            sep = ""))
    ncols <- 9
    if (is.null(x$strata)) {
        if (rmean == "none") 
            end.time <- NA
        else if (is.numeric(rmean)) 
            end.time <- rmean
        else end.time <- max(x$time)
        if (is.matrix(surv)) {
            out <- matrix(0, ncol(surv), ncols)
            for (i in 1:ncol(surv)) {
                if (is.null(x$conf.int)) 
                  out[i, ] <- pfun(x$n, stime, surv[, i], x$n.risk, 
                    x$n.event, NULL, NULL, start.time, end.time)
                else out[i, ] <- pfun(x$n, stime, surv[, i], 
                  x$n.risk, x$n.event, x$lower[, i], x$upper[, 
                    i], start.time, end.time)
            }
            dimnames(out) <- list(dimnames(surv)[[2]], plab)
        }
        else {
            out <- matrix(pfun(x$n, stime, surv, x$n.risk, x$n.event, 
                x$lower, x$upper, start.time, end.time), nrow = 1)
            dimnames(out) <- list(NULL, plab)
        }
    }
    else {
        nstrat <- length(x$strata)
        stemp <- rep(1:nstrat, x$strata)
        last.time <- (rev(x$time))[match(1:nstrat, rev(stemp))]
        if (rmean == "none") 
            end.time <- rep(NA, nstrat)
        else if (is.numeric(rmean)) 
            end.time <- rep(rmean, nstrat)
        else if (rmean == "common") 
            end.time <- rep(median(last.time), nstrat)
        else end.time <- last.time
        if (is.matrix(surv)) {
            ns <- ncol(surv)
            out <- matrix(0, nstrat * ns, ncols)
            if (is.null(dimnames(surv)[[2]])) 
                dimnames(out) <- list(rep(names(x$strata), rep(ns, 
                  nstrat)), plab)
            else {
                cname <- outer(dimnames(surv)[[2]], names(x$strata), 
                  paste, sep = ", ")
                dimnames(out) <- list(c(cname), plab)
            }
            k <- 0
            for (i in 1:nstrat) {
                who <- (stemp == i)
                for (j in 1:ns) {
                  k <- k + 1
                  if (is.null(x$lower)) 
                    out[k, ] <- pfun(x$n[i], stime[who], surv[who, 
                      j], x$n.risk[who], x$n.event[who], NULL, 
                      NULL, start.time, end.time[i])
                  else out[k, ] <- pfun(x$n[i], stime[who], surv[who, 
                    j], x$n.risk[who], x$n.event[who], x$lower[who, 
                    j], x$upper[who, j], start.time, end.time[i])
                }
            }
        }
        else {
            out <- matrix(0, nstrat, ncols)
            dimnames(out) <- list(names(x$strata), plab)
            for (i in 1:nstrat) {
                who <- (stemp == i)
                if (is.null(x$lower)) 
                  out[i, ] <- pfun(x$n[i], stime[who], surv[who], 
                    x$n.risk[who], x$n.event[who], NULL, NULL, 
                    start.time, end.time[i])
                else out[i, ] <- pfun(x$n[i], stime[who], surv[who], 
                  x$n.risk[who], x$n.event[who], x$lower[who], 
                  x$upper[who], start.time, end.time[i])
            }
        }
    }
    if (is.null(x$lower)) 
        out <- out[, 1:7, drop = F]
    if (rmean == "none") 
        out <- out[, -(5:6), drop = F]
    list(matrix = out[, , drop = T], end.time = end.time)
}


