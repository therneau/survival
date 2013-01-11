# Compute the K-M for left/right/interval censored data via Turnbull's
#      slow EM calculation
# x is a factor giving the groups, y is a survival object
survfitTurnbull <- function(x, y, casewt,
                       type=c('kaplan-meier', 'fleming-harrington', 'fh2'),
                       error=c('greenwood', "tsiatis"), se.fit=TRUE,
                       conf.int= .95,
                       conf.type=c('log',  'log-log',  'plain', 'none'),
                       conf.lower=c('usual', 'peto', 'modified'),
		       start.time) {
			     
    type <- match.arg(type)
    method <- match(type, c("kaplan-meier", "fleming-harrington", "fh2"))

    error <- match.arg(error)
    error.int <- match(error, c("greenwood", "tsiatis"))
    conf.type <- match.arg(conf.type)
    conf.lower<- match.arg(conf.lower)
    if (is.logical(conf.int)) {
        # A common error is for users to use "conf.int = FALSE"
        #  it's illegal, but allow it
        if (!conf.int) conf.type <- "none"
        conf.int <- .95
    }

    if (!is.Surv(y)) stop("y must be a Surv object")
    if (!is.factor(x)) stop("x must be a factor")
    xlev <- levels(x)   # Will supply names for the curves
    x <- as.numeric(x)  # keep only the levels

    if (!missing(start.time)) { 
        # The user has requested that survival be "survival given that they
        #  made it to start.time".  We do this by just tossing those who
        #  are known to end before start.time.  Now if one of the times were
        #  interval censored (15,42) and start.time were 20, perhaps it should
        #  be modified too, but we don't.  I really don't know what the 
        #  correct action would be, actually.
        #
	ny <- ncol(y)      
	n.all <- c(table(x))   # remember the original data size
	# remove any obs whose end time is <= start.time
	keep <- (y[,ny-1] >= start.time)
	if (all(keep==FALSE))
		stop(paste("start.time =", start.time,
			   "is greater than all time points."))
	x <- x[keep]
	y <- y[keep,,drop=FALSE]  #make sure y remains a matrix
	casewt <- casewt[keep]
        }
    n.used <- as.vector(table(x))    # This is for the printout
    nstrat <- length(n.used)

    # Make sure that the time variable is not "counting" type,
    #  and convert "left" to "interval" style.
    stype <- attr(y, 'type')
    if (stype=='counting') 
	    stop("survfitTurnbull not appropriate for counting process data")
    if (stype=='interval') status <- y[,3]
    if (stype=='left') status <- ifelse(y[,2]==0,2,1)
    if (stype=='right')status <- y[,2]

    # If any exact times were represented as interval censored, e.g. (x,x)
    #  as the interval for some x, change the code to "uncensored".
    if (any(status==3)) {
        who <- (status==3 & y[,1]==y[,2])
        status[who] <- 1
        }

    #
    # Deal with the issue of "almost tied" times that are equal to
    #  machine precision
    # Convert to the levels of a factor.  We come back to times later
    if (ncol(y)==3) {
        ftime <- factor(y[,1:2])
        y <- matrix(as.numeric(ftime), ncol=2)
    }
    else {
        ftime <- factor(y[,1])
        y <- as.matrix(as.numeric(ftime))
    }

    # the code below actually does the estimate, one curve at a time
    doit <- function(y,status, wt, ...) {
	n <- length(status)
	# Find all of the jump points for the KM in the data set, which are
	#  the exact times, plus any right-followed-by-left pairs.  
        # For this computation, an interval censored observation is considered
        #  to be of the form (a,b], left censored is (-infinity,b] and right
        #  censored is (a, infinity).  If there are two interval censored
        #  obs of (10,20] and (20,40], we do NOT want to create a jtimes entry
        #  at 20.
        #  
	# The algorithm puts a [ at t for each exact, a ( at t for each right
        #  censored, a ] at t for each left censored, and ( and ] at t1/t2 
        #  for each interval censor.  In ties, order the parens as [, ],  (.
        #  Then find pairs of left-followed-immediately-by-right.  The stat2
        #  variable is 0= [ at t, 1= ] at t, 2= ( at t.
        # The variables time2, stat2 are never needed after jtimes has
        #  been created.
	if (any(status==3)) { #interval censored
            stat2 <- c(c(2,0,1,2)[status+1], rep(1, sum(status==3)))
	    time2 <- c(y[,1], y[status==3,2])
	    }
	else {
	    stat2 <-  c(2,0,1)[status+1]
	    time2 <- y[,1]
	    }

	ord <- order(time2, stat2)
	time2 <- time2[ord]
	stat2 <- stat2[ord]
	n2 <- length(time2)
	pairs <- (stat2[-n2]!=1 & stat2[-1]==1)
	jtimes <- c(time2[stat2==0], .5*(time2[-n2] + time2[-1])[pairs])

	#
	# If any of the left censored times are < min(jtime), then treat
	#  them as though they were exact (for now).  The formal MLE
        #  algebra puts all their mass at an arbitray point between the
        #  smallest of such times and -infinity.  
	#
	mintime <- min(jtimes)
	who <- (status==2 & y[,1] < mintime)
	if (any(who)) {
	    status[who] <- 1
	    jtimes <- c(y[who,1], jtimes)
	    }

	# The KM is computed on a fake data set with njump points
	#  standing in for the left and interval censored observations
        # So tempy contains the exact and right censored y data, followed
        #  by the fakes
	jtimes <- sort(unique(jtimes))
	njump <- length(jtimes)
	nreal <- sum(status<2)
	tempx <- factor(rep(1, njump + nreal))  #dummy x var for survfit.km
	tempy <- Surv(c(y[status<2, 1], jtimes),
		      c(status[status<2], rep(1, njump)))

	# wtmat marks, for each left/interval obs, which jump points are in it
        # A column is a "fake" time point, a row is an observation
        # For a left censored obs, we assume that the true event time is 
        #   <= the time recorded, and for an interval one that (a, b] contains
        #   the true event time.  This is motivated by data that would come
        #   from repeated visits, and agrees with Turnbull's paper.
        # If all status <=1, this is the unusual case of left censoring before
        #  some minimal time, in which case I can skip this step.  There are
        #  no interval censored or left censored be "split".
        if (any(status>1)) {
            temp <- matrix(jtimes, nrow=sum(status>1), ncol=njump, byrow=TRUE)
            indx <- (1:n)[status>1] #the subjects of interest
            temp1 <- (temp <= y[indx,1]) # logical matrix for the left censored
            if (any(status>2)) #interval censored
                temp2 <- (temp >  y[indx,1] & temp <= y[indx,2])
            else temp2 <- FALSE & temp1
            temp3 <- rep(status[indx]==2, njump)
            wtmat <- matrix(as.numeric((temp3&temp1) | (!temp3 & temp2)),
                        ncol=njump)
	lwt <- wt[indx]  # the input vector of case weights, for these
        }
        else {
            wtmat <- matrix(rep(1, length(jtimes)), nrow=1)
            lwt <- 1
        }

	eps <- 1

	# The initial "starter" KM is proportional to the number of intervals
        #  that overlap each time point
        temp <- apply(wtmat, 2, sum)
	tfit  <- list(time=jtimes, surv= 1- cumsum(temp)/sum(temp))
	old <- tfit$surv

        iter <- 0
        aitken1 <- jump1 <- 0 #dummy values for lagging
	while (eps > .00005) {
            iter <- iter +1
	    # partition each left/interval person out over the jumps
	    jumps <- -diff(c(1, tfit$surv[match(jtimes, tfit$time)])) #KM jumps
            if (TRUE) { # add Aitken acceleration to speed things up
                # Given 3 points on a sequence, it guesses ahead.  So we use
                #  a set of 3 to guess ahead, generate 3 more regular EM,
                #  guess ahead, 3 regular EM, etc. 
                # Actually, we go every 5th below instead of every 3rd, to
                #  give the EM a chance to restabilize the relative
                #  sizes of elements of "jump".  We also only allow it to
                #  stop when comparing two "real EM" iterations.
                aitken2 <- aitken1      
                aitken1 <- jumps - jump1
                jsave <- jumps
                if (iter%%5 ==0) {
                    oldlik <- sum(log(wtmat %*% jumps))
                    jumps <- jump2 - (aitken2)^2/(aitken1 - aitken2)
                    bad <- (jumps<=0 | jumps >=1 | is.na(jumps))
                    jumps[bad] <- jsave[bad]  #failsafe
                    newlik <- sum(log(wtmat %*% jumps))
                    if (newlik < oldlik) jumps <- jsave # aitkin didn't work!
                    }
                jump2   <- jump1  # jumps, lagged by 2 iterations
                jump1   <- jsave  # jumps, lagged by 1 iteration
                }

	    wt2 <- wtmat %*% diag(jumps)
	    wt2 <- (lwt/(apply(wt2,1,sum))) * wt2 
	    wt2 <- apply(wt2, 2, sum)
	    tfit <- survfitKM(tempx, tempy, casewt=c(wt[status<2], wt2), ...)

	    if (FALSE) {
                # these lines are in for debugging: change the above to 
                #  " if (TRUE)" to turn on the printing
                cat("\n Iteration = ", iter, "\n")
		cat("survival=",
		    format(round(tfit$surv[tfit$n.event>0],3)),  "\n")
		cat(" weights=", format(round(wt2,3)), "\n")
		}
            stemp <- tfit$surv[match(jtimes, tfit$time)] 
            if (iter%%5<2) eps <- 1  #only check eps for a pair of EM iters
            else eps <- max(abs(old-stemp))
	    old <- stemp
	    }	

	#
	# Now, fix up the "cheating" I did for any left censoreds which were
	#  less than the smallest jump time
	who <- (tfit$time < mintime & tfit$n.event >0)
	if (any(who)) {
	    indx <- match(mintime, tfit$time)  # first "real" time
#	    tfit$surv[who] <- tfit$surv[indx]
	    tfit$n.event[who] <- 0
#	    if (!is.null(tfit$std.err)) {
#		tfit$std.err[who] <- tfit$std.err[indx]
#		tfit$lower[who]   <- tfit$lower[indx]
#		tfit$upper[who]   <- tfit$upper[indx]
#		}
	    }
	tfit
	}		
    #
    # Now to work, one curve at a time
    #
    time   <- vector('list', nstrat)
    n.risk <- vector('list', nstrat)
    surv   <- vector('list', nstrat)
    n.cens <- vector('list', nstrat)
    n.event<- vector('list', nstrat)
    fnumeric <- as.numeric(levels(ftime))

    uniquex <- sort(unique(x))
    for (i in 1:nstrat) {
	who <- (x== uniquex[i])
	tfit <- doit(y[who,,drop=FALSE], status[who], casewt[who])
	time[[i]]   <- approx(1:length(fnumeric), fnumeric, tfit$time)$y
	n.risk[[i]] <- tfit$n.risk
	surv[[i]]   <- tfit$surv
	n.cens[[i]] <- tfit$n.cens
	n.event[[i]]<- tfit$n.event
	if (i==1) {
	    if (!is.null(tfit$std.err)) {
		std.err <- vector('list', nstrat)
		conf.lower <- vector('list', nstrat)
		conf.upper <- vector('list', nstrat)
		se.fit <- TRUE
		}
	    else se <- FALSE
	    }
	if (se.fit) {
	    std.err[[i]]    <- tfit$std.err
	    conf.lower[[i]] <- tfit$lower
	    conf.upper[[i]] <- tfit$upper
	    }
	}

    temp <- list(n=n.used,
		 time = unlist(time),
		 n.risk = unlist(n.risk),
		 n.event= unlist(n.event),
		 n.censor = unlist(n.cens),
		 surv = unlist(surv),
		 type='right')
    
    if (nstrat >1) {
        strata <- unlist(lapply(time, length))
	names(strata) <- xlev[sort(unique(x))]
	temp$strata <- strata
	}

    if (se.fit) {
	temp$std.err <- unlist(std.err)
	temp$lower <- unlist(conf.lower)
	temp$upper <- unlist(conf.upper)
	temp$conf.type <- tfit$conf.type
	temp$conf.int  <- tfit$conf.int
	}
    temp
    }
