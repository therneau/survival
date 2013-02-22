#  A version that does more work in S, less in C
survfitKM <- function(x, y, casewt=rep(1,length(x)),
		       type=c('kaplan-meier', 'fleming-harrington', 'fh2'),
		       error=c('greenwood', "tsiatis"), se.fit=TRUE,
		       conf.int= .95,
		       conf.type=c('log',  'log-log',  'plain', 'none'),
		       conf.lower=c('usual', 'peto', 'modified'),
		       start.time, new.time) {
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
    if (attr(y, 'type') != 'right' && attr(y, 'type') != 'counting')
	    stop("Can only handle right censored or counting data")
    ny <- ncol(y)       # Will be 2 for right censored, 3 for counting
    xlev <- levels(x)   # Will supply names for the curves
    x <- as.numeric(x)  # keep only the levels

    #  Allow "new.time" as a synonym for start.time
    if (missing(start.time) && !missing(new.time)) start.time <- new.time
    if (!missing(start.time)) { 
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

    #
    # Each of the necessary output objects is originally a list with one
    #  element per strata.  This doesn't use up extra S memory, the number
    #  of curves is usually small enough that the "for" loop is no great
    #  cost, and it's easier to see what's going on than C code.
    # If tapply gets fast again, this code will also be fast.
    # The remaining C-code is simple, but hard to do cleanly in S
    #  Let nrisk=A and nevent=B.  The terms in returned sum1 and sum2 are
    #     If ndead= 0, sum1=sum2=1 (avoids a 0/0 in the S code)
    #     If ndead= 1, sum1=1/A and sum2= 1/A*A
    #     If ndead= 2, sum1= (1/2)[ 1/A + 1/(A- B/2)]
    #                  sum2= (1/2)[ 1/A^2 + (1/(A-B/2))^2]
    #     If ndead =3, sum1 = (1/3)[1/A + 1/(A-B/3) + 1/(A -2B/3)]
    #     If ndead =4, sum1 = (1/4)[1/A + 1/(A-B/4) + 1/(A -2B/4) + 1/(A-3B/4)]
    #		 and etc.
    time   <- vector('list', nstrat)
    n.risk <- vector('list', nstrat)
    surv   <- vector('list', nstrat)
    n.cens <- vector('list', nstrat)
    n.event<- vector('list', nstrat)
    strata <- integer(nstrat)
    if (se.fit) varhaz <- vector('list', nstrat)
    if (ny==3)  n.enter <- vector('list', nstrat)

    uniquex <- sort(unique(x))
    for (i in 1:nstrat) {
	who <- (x== uniquex[i])
	if (ny==2) {
	    # the "factor" + levels ensures 2 columns in temp even if all 
	    #  are dead or  all are alive
	    temp <- tapply(casewt[who], 
			   list(factor(y[who,1]), 
                                factor(y[who,2], levels=0:1)), sum)
	    temp <- ifelse(is.na(temp), 0, temp)

            # The two lines below do not always give the same answer
            #  When two times differ by only the machine precision, unique
            #  will give more values, and then "time" will be a different
            #  length than the other components
	    time[[i]] <- as.numeric(dimnames(temp)[[1]])
#	    time[[i]] <- sort(unique(y[who,1]))   #same answer as line above
	    ntemp  <- (dim(temp))[1]
	    nevent <- as.vector(temp[,2])
	    ncens  <- as.vector(temp[,1])
	    nrisk  <- rev(cumsum(rev(temp %*% c(1,1))))
	    ndead  <- as.vector(table(y[who,1], factor(y[who,2], 
                                                       levels=0:1)) [,2])
	    }
	else {
	    #    The counting process case
	    # We have to be a bit more clever here -- if I did a table of the
	    #  start times and a separate one of the stop times, they wouldn't
	    #  necessarily match.  So do it all at once with a fake 'status'
	    #  variable which is ==2 for start times and = status for stop
	    # Also, the number of rows in output (ntemp) may be bigger than
	    #  the number of rows of input (n).
	    n <- sum(who)
	    temp <- factor(c(rep(2,n),y[who,3]), levels=0:3)
	    temp <- tapply(rep(casewt[who],2), list(factor(y[who,1:2]), temp),
                           sum)
	    temp <- ifelse(is.na(temp), 0, temp)
	    time[[i]] <- as.numeric(dimnames(temp)[[1]])
	    ntemp  <- (dim(temp))[1]
	    nevent <- as.vector(temp[,2])
	    ncens  <- as.vector(temp[,1])
	    nenter <- as.vector(temp[,3])
	    nrisk  <- cumsum(nenter - (nevent + ncens))
	    nrisk  <- c(0, nrisk[-ntemp])  #risk counts change at time t+0
	    n.enter[[i]] <- nenter
	    # again, a fake status to make sure that all the times appear
	    ndead <- as.vector(table(y[who,1:2], 
				 factor(c(rep(0,n),y[who,3]), levels=0:2))[,2])
	    }
	
	strata[i] <- ntemp
	trisk <- ifelse(nrisk==0, 1, nrisk) #avoid 0/0 cases

	if (method==1) surv[[i]] <- cumprod((trisk-nevent)/trisk)
	if (method==2) {
	    hazard <- nevent/trisk  #Nelson's hazard estimate
	    surv[[i]] <- exp(- cumsum(hazard))
	    }
	if (method==3) {
	    tsum  <- .C(Csurvfit4, as.integer(length(ncens)),
			            as.integer(ndead),
			            sum1 = as.double(nrisk),
			            sum2 = as.double(nevent))
	    hazard <- nevent *tsum$sum1
	    surv[[i]] <- exp(-cumsum(hazard))
	    }

	if (se.fit) {
	    if (error.int==1) # Greenwood
		    varhaz[[i]] <- cumsum(nevent/(trisk*(trisk-nevent)))
	    else {
		if (method <3) varhaz[[i]] <- cumsum(nevent/(trisk^2))
		else 	       varhaz[[i]] <- cumsum(nevent* tsum$sum2)
		}
	    }
	n.event[[i]] <- nevent
	n.cens[[i]]  <- ncens
	n.risk[[i]]  <- nrisk
	}

    if (ny==2) {
	temp <- list(n=n.used,
		     time = unlist(time),
		     n.risk = unlist(n.risk),
		     n.event= unlist(n.event),
		     n.censor = unlist(n.cens),
		     surv = unlist(surv),
		     type='right')
	}
    else {
	temp <- list(n=n.used,
		     time = unlist(time),
		     n.risk = unlist(n.risk),
		     n.event= unlist(n.event),
		     n.censor = unlist(n.cens),
		     n.enter  = unlist(n.enter),
		     surv = unlist(surv),
		     type='counting')
	}

    if (nstrat >1) {
	names(strata) <- xlev[sort(unique(x))]
	temp$strata <- strata
	}

    if (!missing(start.time)) {
	    temp$start.time <- start.time # user defined time to start
	    temp$n.all <- n.all
	    }

    if (se.fit) {
	std.err <- sqrt(unlist(varhaz))
	temp$std.err <- std.err
	#
	# n.lag = the # at risk the last time there was an event (or
	#   the first time of a strata)
	#
	events <- temp$n.event >0
	if (nstrat==1) events[1] <- TRUE
	else           events[1 + cumsum(c(0, strata[-nstrat]))] <- TRUE
	zz <- 1:length(events)
	n.lag <- rep(temp$n.risk[events], diff(c(zz[events], 1+max(zz))))
	std.low <- switch(conf.lower,
			  'usual' = std.err,
			  'peto' = sqrt((1-temp$surv)/ temp$n.risk),
			  'modified' = std.err * sqrt(n.lag/temp$n.risk))
	zval <- qnorm(1- (1-conf.int)/2, 0,1)

	if (conf.type=='plain') {
	    temp1 <- temp$surv + zval* std.err * temp$surv
	    temp2 <- temp$surv - zval* std.low * temp$surv
	    temp <- c(temp, list(upper=pmin(temp1,1), lower=pmax(temp2,0),
				 conf.type='plain', conf.int=conf.int))
	    }

	if (conf.type=='log') {
	    #avoid some "log(0)" messages
	    xx <- ifelse(temp$surv==0,1,temp$surv)  

	    temp1 <- ifelse(temp$surv==0, NA, exp(log(xx) + zval* std.err))
	    temp2 <- ifelse(temp$surv==0, NA, exp(log(xx) - zval* std.low))
	    temp <- c(temp, list(upper=pmin(temp1,1), lower=temp2,
				 conf.type='log', conf.int=conf.int))
	    }

	if (conf.type=='log-log') {
	    who <- (temp$surv==0 | temp$surv==1) #special cases
	    temp3 <- ifelse(temp$surv==0, NA, 1)
	    xx <- ifelse(who, .1,temp$surv)  #avoid some "log(0)" messages
	    temp1 <- exp(-exp(log(-log(xx)) + zval*std.err/log(xx)))
	    temp1 <- ifelse(who, temp3, temp1)
	    temp2 <- exp(-exp(log(-log(xx)) - zval*std.low/log(xx)))
	    temp2 <- ifelse(who, temp3, temp2)
	    temp <- c(temp, list(upper=temp1, lower=temp2,
				 conf.type='log-log', conf.int=conf.int))
	    }
        }

    temp
    }
