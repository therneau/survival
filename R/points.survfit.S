# $Id: points.survfit.S 11390 2010-02-24 19:22:25Z therneau $
points.survfit <- function(x, xscale, xmax, fun, ...) {
    ssurv <- x$surv
    stime <- x$time

     if (inherits(x, "survfitms")) {
        x$surv <- 1- x$prev
        if (is.matrix(x$surv)) dimnames(x$surv) <- list(NULL, x$states)
        if (!is.null(x$lower)) {
            x$lower <- 1- x$lower
            x$upper <- 1- x$upper
        }
        if (missing(fun)) fun <- "event"
    }

   if (is.null(x$strata)) {
	nstrat <- 1
	stemp <- rep(1, length(x$time))
	}
    else {
	nstrat <- length(x$strata)
	stemp <- rep(1:nstrat, x$strata)
	}
    if (!missing(xmax) && any(stime>xmax)) {
	stime[stime>xmax] <- NA # simple way to stop it from plotting
	}

    stime <- stime/xscale

    if (!missing(fun)) {
        if (is.character(fun)) {
            tfun <- switch(fun,
                            'log' = function(x) x,
                            'event'=function(x) 1-x,
                            'cumhaz'=function(x) -log(x),
                            'cloglog'=function(x) log(-log(x)),
                            'pct' = function(x) x*100,
                            'logpct'= function(x) 100*x,
                            stop("Unrecognized function argument")
                            )
            }
        else if (is.function(fun)) tfun <- fun
        else stop("Invalid 'fun' argument")

        ssurv <- tfun(ssurv)
        }


    if (!is.matrix(ssurv))
	    points(stime, ssurv, ...)
    else
	    matpoints(stime, ssurv, ...)
    }
