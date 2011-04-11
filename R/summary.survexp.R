#
# Almost identical to summary.survfit.  The big differences
#  are no call to the survmean function (irrelevant), and
#  there is no censoring, extend, or rmean argument.
# Because survexp objects do not contain n.event or n.censor, 
#  subsetting is easier.

summary.survexp <- function(object, times, scale=1, ...) {
    fit <- object
    if (!inherits(fit, 'survexp'))
	    stop("Invalid data")

    # The fit$surv object is sometimes a vector and sometimes a matrix.
    #  Make a copy of it that is always a matrix, to simplify the number of
    #  cases for our subscripting work below.  At the end of the routine
    #  we'll turn it back into a vector if needed.  Similar treatment is
    #  given to the standard error and confidence limits.
    surv <- as.matrix(fit$surv)
    if (is.null(fit$strata)) {
	nstrat <- 1
	stemp <- rep(1, nrow(surv))
        strata.names <- ""
	}
    else   {
	nstrat <- length(fit$strata)
	stemp <- rep(1:nstrat, fit$strata)
        strata.names <- names(fit$strata)
	}

    if (is.null(fit$std.err)) std.err <- NULL
    else 		      std.err <- fit$std.err * surv

    if (!is.null(fit$lower)) {
	lower <- as.matrix(fit$lower)
	upper <- as.matrix(fit$upper)
        }

    if (missing(times)) {
        times <- fit$time
        n.risk<- fit$n.risk
        strata <- factor(stemp, labels=strata.names)
    }
    else {  
	#this case is harder, since it involves "in between" points
	times <- sort(times)   #just in case the user didn't

	# The basic idea is to process the curves one at a time,
	#   adding the results for that curve onto a list, so the
	#   survival surv[[1], surv[[2]], etc.
	# For the survival, stderr, and confidence limits it suffices
	#   to create a single list 'indx1' containing a subscripting vector
	indx1 <- n.risk <- newtimes <- vector('list', nstrat)
	n <- length(stemp)
	for (i in 1:nstrat) {
	    who <- (1:n)[stemp==i]  # the rows of the object for this strata
	    stime <- fit$time[who]

	    # First, toss any printing times that are outside our range
	    mintime <- min(stime, 0)
	    ptimes <- times[times >= mintime]
            maxtime <- max(stime)
            ptimes <- ptimes[ptimes <= maxtime]

	    newtimes[[i]] <- ptimes

	    # If we tack a -1 onto the front of the vector of survival
	    #  times, then indx1 is the subscript for that vector
	    #  corresponding to the list of "ptimes".  If the input
	    #  data had stime=c(10,20) and ptimes was c(5,10,15,20),
	    #  the result would be 1,2,2,3.
	    # For n.risk we want a slightly different index: 2,2,3,3.
	    #  "In between" times point to the next higher index for n.risk,
	    #  but the next lower one for survival. (Survival drops at time t,
	    #  the n.risk immediately afterwords at time t+0: you were at
	    #  risk just before you die, but not a moment after). The
	    #  extra point needs to be added at the end.
	    #
	    ntime <- length(stime)  #number of points
	    temp1 <- approx(c(mintime-1, stime), 0:ntime, xout=ptimes,
			    method='constant', f=0, rule=2)$y
	    indx1[[i]] <- ifelse(temp1==0, 1, 1+ who[pmax(1,temp1)])
            # Why not just "who[temp1]" instead of who[pmax(1,temp1)] in the
            #  line just above?  When temp1 has zeros, the first expression
            #  gives a vector that is shorter than temp1, and the ifelse
            #  doesn't work right due to mismatched lengths.  

	    # Compute the number at risk.  If stime = 1,10, 20 and ptime=3,10,
	    #   12, then temp1 = 2,2,3: the nrisk looking ahead
	    # approx() doesn't work if stime is of length 1
	    if (ntime ==1) temp1 <- rep(1, length(ptimes))
	    else temp1 <- approx(stime, 1:ntime, xout=ptimes,
			    method='constant', f=1, rule=2)$y
	    n.risk[[i]] <- ifelse(ptimes>max(stime), 
                                  fit$n.risk[length(fit$n.risk)],
				  fit$n.risk[who[temp1]])
	    }

	# Now create the output list
	times  <- unlist(newtimes)
	n.risk <-  unlist(n.risk)

	indx1 <- unlist(indx1)
	surv <- (rbind(1.,surv))[indx1,,drop=FALSE]
	if (!is.null(std.err)) std.err <- rbind(0.,std.err)[indx1,,drop=FALSE]
	if (!is.null(fit$lower)) {
	    lower <- rbind(1.,lower)[indx1,,drop=FALSE]
	    upper <- rbind(1.,upper)[indx1,,drop=FALSE]
	    }
	if (!is.null(fit$strata)) {
	    scount <- unlist(lapply(newtimes, length))
	    strata <- factor(rep(1:nstrat, scount), labels=names(fit$strata))
	    }
	}

    #
    # Final part of the routine: paste the material together into
    #  the correct output structure
    #
    temp <- list(surv=surv, time=times/scale, n.risk=n.risk)

    if (ncol(surv)==1) {
	# Make surve & etc vectors again
	temp$surv <- drop(temp$surv)
	if (!is.null(std.err)) temp$std.err <- drop(std.err)
	if (!is.null(fit$lower)) {
	    temp$lower <- drop(lower)
	    temp$upper <- drop(upper)
	    }
	}
    else {
	if (!is.null(std.err)) temp$std.err <- std.err
	if (!is.null(fit$lower)) {
	    temp$lower <- lower
	    temp$upper <- upper
	    }
	}

    if (!is.null(fit$strata)) {
	temp$strata <- strata
	}
    
    temp$method <- fit$method
    temp$call <- fit$call
    if (!is.null(fit$na.action)) temp$na.action <- fit$na.action
  
    class(temp) <- "summary.survexp"
    temp
    }









