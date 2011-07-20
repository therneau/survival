#
# Compute the current incidence curve for a data set
#   called from survfit.formula
survfitCI <- function(X, Y, weights, 
                       type=c('kaplan-meier', 'fleming-harrington', 'fh2'),
                       error=c('greenwood', "tsiatis"), se.fit=TRUE,
                       conf.int= .95,
                       conf.type=c('log',  'log-log',  'plain', 'none'),
                       conf.lower=c('usual', 'peto', 'modified'),
                       etype, id) {
    # Turn off standard errors, as there is a major bug 10/2010
    se.fit <- FALSE

    method <- match.arg(type)
    error <- match.arg(error)
    error.int <- match(error, c("greenwood", "tsiatis"))
    conf.type <- match.arg(conf.type)
    conf.lower<- match.arg(conf.lower)

    type <- attr(Y, "type")
    if (type!='right' && type!='counting')
	stop(paste("Cumulative incidence computation doesn't support \"", type,
			  "\" survival data", sep=''))

    # Now to work
    #       Surv(time, status) ~ sex, etype= event.type
    # is the canonical call.  This will produce a matrix of survivals,
    # stratified by sex.  Optionally, there can be an id statement
    # to indicate a multiple transitions/subject data set.
    #
    #   If status==1, then a transition at time "time" is assumed (or at stop
    # time for (start, stop] data.  If it is 0, then the subject is
    # censored.  If censored, the value of the event type variable
    # is completely ignored.  However, don't use NA in that case, since
    # the na.action will remove that line of data!  Subjects are assumed
    # to start in a "null" state, which is not tabulated for survival.
    # To change this behavior, give subjects a transition at time 0.
    #   
    # If there is not a id term, then the result can be computed
    #  by post-processing the KM.  If there is a id term, then we
    #  must explicitly use the redistribut-to-the-right algorithm to 
    #  compute the estimate
    #
    # Due to roundoff, unique(time) and table(time) may not be the same.
    #   Dodge the issue by turning it into a factor.
    ftime <- factor(Y[,-ncol(Y)])
    Y[, -ncol(Y)] <- as.numeric(ftime)
    time <- Y[,ncol(Y)-1]  #unique death times
    status <- Y[,ncol(Y)]
    ncurve <- length(levels(X))
    n <- length(time)

    # make sure that the etypes for status==0 are not tabulated
    #  in our output
    if (is.factor(etype)) {
        temp <- levels(etype[status !=0, drop=TRUE])
        etype <- factor(ifelse(status==0, NA, etype))
        levels(etype) <- temp
        }
    else etype <- factor(ifelse(status==0, NA, etype))
    
    if (missing(id) || length(id)==0 || length(unique(id))==n) {
        # this case implicitly assumes one transition/subject
        # use survfit.km to get the # at risk and etc
	kfit <- survfitKM(X, Y, se.fit=FALSE)
        
        # partition out the event type
	if (is.null(kfit$strata)) {
	    temp <- table(time, etype, factor(status, levels=0:1), 
                          exclude=NA)[,,2]
	    temp <- temp / ifelse(kfit$n.event==0, 1, kfit$n.event)  # percents
	    jumps <- -diff(c(1, kfit$surv))
	    newjump <- jumps * as.matrix(temp) #as.matrix in case only 1 event
                                               # type
	    newsurv <- apply(newjump, 2, function(x) 1 - cumsum(x))
	    }
	else {  # do it one curve at a time
	    netype <- length(levels(etype))
	    newsurv <- matrix(1.0, nrow=length(kfit$time), ncol=netype)
	    newjump <- 0*newsurv
	    stemp <- names(kfit$strata)
	    istart <- cumsum(c(1,kfit$strata)) #starting points
	    for (i in 1:length(kfit$strata)) {
		# j will index the rows in the survival (kfit) object
		# who will index into the original data
		j <- seq(from=istart[i], length=kfit$strata[i])
		who <- (X==stemp[i] & status==1)
		if (any(who)) { # if any deaths
		    temp <- table(factor(time[who], levels=kfit$time[j]),
				  etype[who], exclude=NA)
		    tevent <- kfit$n.event[j]
		    temp <- temp / ifelse(tevent<1, 1, tevent)  # percents
		    jumps <- -diff(c(1, kfit$surv[j]))
		    jumps <- jumps * temp
		    newsurv[j,] <- apply(jumps, 2, function(x) 1 - cumsum(x))
		    newjump[j,] <- jumps  # save for se computations
		    }
		}
	    }
	dimnames(newsurv) <- list(NULL, levels(etype))
	kfit$surv <- newsurv
		
	# compute the Greenwood error
	if (se.fit) {
	    greenwood <- newjump /(kfit$n.risk*(kfit$n.risk-newjump))
	    kfit$std.err <- sqrt(apply(greenwood, 2, cumsum))
	    }
        }
    
    else {
        # Have to do real work in this case
 	if (ncol(Y) !=2) {
	    stop("(start, stop] data for repeated events not yet allowed")
	    }
	time   <- Y[,1]
	status <- Y[,2]

	# ensure consistency of the weights
        idlist <- sort(unique(id))
	firstone <- match(idlist, id) # first obs for each subject
        # if the ids were Smith, Jones, Doe, Smith, Doe, then indx
        #   would be 1,2,3,1,3, ...  pointers to the first obs of each id
        indx <- firstone[match(id, idlist)] 

	if (any(weights[indx] != weights))
		stop("All case weights for a given subject must agree")
	if (ncurve >1) {
	    # Ensure that any given subject is in only one curve
            if (any(X[indx] != X))
		stop("No subject can be in more than 1 curve")
	    }

	# Do some preprocessing, to make things nice integers (1,2, etc)
	etype2 <- as.numeric(etype)      # numeric etypes
        netype <- length(levels(etype))

	# This function creates a single survival curve
	docurve <- function(time, status, etype, weights, id, netype) {
	    n <- length(time)
	    timelist <- sort(unique(time))   # unique event/censor times
	    ntime <- length(timelist)
	    surv   <- matrix(0., nrow=ntime, ncol=netype)
	    varhaz <- surv
	    nrisk <- nevent <- double(ntime)

	    # create a unique weight for each subject
	    # cetype is the current etype for each subject, initially 0 or
	    #  "no etype at all".  The factor causes all survival matrices
	    #  to be the same size.
	    idlist <- unique(id)
	    wtindex <- match(id, idlist)
	    wt2 <- weights[match(idlist, id)]
	    cetype <- factor(rep(0, length(idlist)), levels=0:netype)
            islast <- c(id[-1]!= id[-n], TRUE)  # mark the last obs of each id

            # For computing "number at risk", I will need the original
            #  weights (censorings not redistributed)
            wt.orig <- wt2

	    # Walk down the timelist, one at a time
	    for (i in 1:ntime) {
		xtime <- timelist[i]
                at.risk <- unique(wtindex[time>= xtime])
		nrisk[i] <- sum(wt.orig[at.risk]) #still at risk for transition
		event <- (time==xtime & status >0)
		nevent[i] <- sum(wt.orig[wtindex[event]])
            
		# reassign the etype of those with an event and 
		#  recompute survival
		if (any(event)) {
		    cetype[wtindex[event]] <- etype[event]
		    temp <- tapply(wt2, cetype, sum)/sum(wt2)
		    surv[i,] <- ifelse(is.na(temp), 0, temp)[-1]
		    varhaz[i,] <- (surv[i,]/ (1-surv[i,])) * sum(wt2^2)/
			(sum(wt2))^2
		    }
		else {
		    if (i>1) {
			surv[i,] <- surv[i-1,]
			varhaz[i,] <- varhaz[i-1,]
			}
		    }
		# Now, redistribute the weights if those who are at the end
                # of their last interval.  Note that a death will never get
                # redistributed, because there is no one to recieve it: that
                # would require the data to have an obs with futime after death
		redo <- which(time==xtime & islast) 
                for (j in redo) {
                    # list of others in the same etype, at this time, with
                    #   a positive weight
                    k <- wtindex[j]   # the position of this id in wt2 vector
                    child <- (cetype==cetype[k] & wt2 >0)
                    otherwt <- sum(wt2[child]) - wt2[k]
                    # if no children, the weight is not redistributed
                    if (otherwt >0){
                        wt2[child] <- wt2[child]*(1 + wt2[k]/otherwt)
                        wt2[k] <- 0
                        }
		    }
		}
	    list(time=timelist, surv= 1.0 - surv, n.risk=nrisk,
		 n.event=nevent, std.err= sqrt(varhaz))
	    }
	
        # Ok, now we're actually ready to do the work.  Copy the
	#  list/unlist trick of survfit.km, although the matrices have
	#  to be transposed in order to make it work.
	#
	if (ncurve==1) {
	    kfit <- docurve(time, status, etype2, weights, id, netype)
	    kfit$type <- 'right'
	    kfit$call <- 'call'
	    kfit$n <- length(Y)
            
	    }
	else {
	    timex  <- vector('list', ncurve)
	    n.risk <- vector('list', ncurve)
	    surv   <- vector('list', ncurve)
	    n.event<- vector('list', ncurve)
	    strata <- integer(ncurve)
	    if (se.fit) varhaz <- vector('list', ncurve)
	    uniquex <- levels(X)
	    for (i in 1:ncurve) {
		who <- (X== uniquex[i])
		temp <- docurve(time[who], status[who], etype2[who],
				weights[who], id[who], netype)
		timex[[i]]   <- temp$time
		n.risk[[i]]  <- temp$n.risk
		n.event[[i]] <- temp$n.event
		surv[[i]]    <- c(t(temp$surv))
		strata[i]    <- length(temp$time)
		if (se.fit) varhaz[[i]] <- c(t(temp$std.err))
		}
	    names(strata) <- uniquex
	    kfit <- list(n=table(X),
			 time= unlist(timex),
			 n.risk= unlist(n.risk), 
			 n.event=unlist(n.event),
			 surv = t(matrix(unlist(surv), nrow=netype)),
			 type='right', call=call)
	    if (se.fit) kfit$std.err <- t(matrix(unlist(varhaz),nrow=netype))
	    }
	dimnames(kfit$surv) <- list(NULL, levels(as.factor(etype)))
        }	

    #	
    # Last bit: add in the confidence bands (also stolen from survfit.km)
    #
    if (se.fit) {
	std.err <- kfit$std.err
	#
	# n.lag = the # at risk the last time there was an event (or
	#   the first time of a strata)
	#
	events <- kfit$n.event >0
	if (ncurve==1) events[1] <- TRUE
	else           events[1 + cumsum(c(0, kfit$strata[-ncurve]))] <- TRUE
	zz <- 1:length(events)
	n.lag <- rep(kfit$n.risk[events], diff(c(zz[events], 1+max(zz))))
	std.low <- switch(conf.lower,
			  'usual' = std.err,
			  'peto' = sqrt((1-kfit$surv)/ kfit$n.risk),
			  'modified' = std.err * sqrt(n.lag/kfit$n.risk))
	zval <- qnorm(1- (1-conf.int)/2, 0,1)

	if (conf.type=='plain') {
	    temp1 <- kfit$surv + zval* std.err * kfit$surv
	    temp2 <- kfit$surv - zval* std.low * kfit$surv
	    kfit <- c(kfit, list(upper=pmin(temp1,1), lower=pmax(temp2,0),
				 conf.type='plain', conf.int=conf.int))
	    }

	if (conf.type=='log') {
	    #avoid some "log(0)" messages
	    xx <- ifelse(kfit$surv==0,1,kfit$surv)  

	    temp1 <- ifelse(kfit$surv==0, NA, exp(log(xx) + zval* std.err))
	    temp2 <- ifelse(kfit$surv==0, NA, exp(log(xx) - zval* std.low))
	    kfit <- c(kfit, list(upper=pmin(temp1,1), lower=temp2,
				 conf.type='log', conf.int=conf.int))
	    }

	if (conf.type=='log-log') {
	    who <- (kfit$surv==0 | kfit$surv==1) #special cases
	    temp3 <- ifelse(kfit$surv==0, NA, 1)
	    xx <- ifelse(who, .1,kfit$surv)  #avoid some "log(0)" messages
	    temp1 <- exp(-exp(log(-log(xx)) + zval*std.err/log(xx)))
	    temp1 <- ifelse(who, temp3, temp1)
	    temp2 <- exp(-exp(log(-log(xx)) - zval*std.low/log(xx)))
	    temp2 <- ifelse(who, temp3, temp2)
	    kfit <- c(kfit, list(upper=temp1, lower=temp2,
				 conf.type='log-log', conf.int=conf.int))
	    }
        }

    # Now replace the integer times with the original time scale
    temp <- as.numeric(levels(ftime))
    kfit$time <- temp[kfit$time]
    kfit
    }

