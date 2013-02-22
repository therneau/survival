coxph.detail <-  function(object, riskmat=FALSE) {
    method <- object$method
    if (method!='breslow' && method!='efron')
	stop(paste("Detailed output is not available for the", method,
			"method"))
    n <- length(object$residuals)
    rr <- object$residuals
    weights <- object$weights        #always present if there are weights
    x <- object[['x']]
    y <- object$y
    strat <- object$strata
    Terms <- object$terms
    if (!inherits(Terms, 'terms'))
	    stop("invalid terms component of object")
    strats <- attr(Terms, "specials")$strata

    if (is.null(y)  ||  is.null(x)) {
        mf <- model.frame(object)
        y <- model.response(mf)
        x <- model.matrix(object, data=mf)
        if (length(strats)) {
            stemp <- untangle.specials(object$terms, 'strata', 1)
            if (length(stemp$vars)==1) strat <- mf[[stemp$vars]]
            else strat  <- strata(mf[,stemp$vars], shortlabel=TRUE)
            }
	}

    nvar <- ncol(x)
    if (ncol(y)==2) {
	mintime <- min(y[,1])
	if (mintime < 0) y <- cbind( 2*mintime -1, y)
	else 	y <- cbind(-1,y)
	}
    if (is.null(strat)) {
	ord <- order(y[,2], -y[,3])
	newstrat <- rep(0,n)
	}
    else {
	ord <- order(strat, y[,2], -y[,3])
	newstrat <- c(diff(as.numeric(strat[ord]))!=0 ,1)
	}
    newstrat[n] <- 1

    # sort the data
    x <- x[ord,]
    y <- y[ord,]
    storage.mode(y) <- 'double'
    score <- exp(object$linear.predictors)[ord]
    if (is.null(weights)) weights <- rep(1,n)
    else                  weights <- weights[ord]

    ndeath <- sum(y[,3])
    if (riskmat) {
	rmat <- integer(ndeath*n)
	}
    else rmat <- as.integer(1)
    
    
    ff <- .C(Ccoxdetail, as.integer(n),
			  as.integer(nvar),
			  ndeath= as.integer(ndeath),
			  y = y,
			  as.double(x),
			  index = as.integer(newstrat),
			  event2 =as.double(score),
			  weights = as.double(weights),
			  means= c(method=='efron', double(ndeath*nvar-1)),
			  u = double(ndeath*nvar),
			  i = double(ndeath*nvar*nvar),
	                  rmat = rmat,
	                  nrisk2 = double(ndeath),
			  double(nvar*(3 + 2*nvar)))
     keep <- 1:ff$ndeath
    vname<- dimnames(x)[[2]]
    time <- y[ff$index[keep],2]
    names(time) <- NULL
    means<- (matrix(ff$means,ndeath, nvar))[keep,]
    score<-  matrix(ff$u, ndeath, nvar)[keep,]
    var <- array(ff$i, c(nvar, nvar, ndeath))[,,keep]
    if (riskmat) {
	rmat <- matrix(0, n, ff$ndeath)
	rmat[ord,] <- ff$rmat[1:(n*ff$ndeath)]  # in the order of orig data
	dimnames(rmat) <- list(NULL, time)
	}

    if (nvar>1) {
	dimnames(means) <- list(time, vname)
	dimnames(score) <- list(time, vname)
	dimnames(var) <- list(vname, vname, time)
	}
    else {
	names(means) <- time
	names(score) <- time
	names(var) <- time
	}

    dimnames(ff$y) <- NULL
    temp <- list(time = time, means=means, nevent=ff$y[keep,1],
	 nrisk = ff$y[keep,2], hazard= ff$y[keep,3], score= score,  imat=var,
	 varhaz=ff$weights[keep], y=y, x=x)
    if (length(strats)) temp$strata <- table((strat[ord])[ff$index[keep]])
    if (riskmat) temp$riskmat <- rmat
    if (!all(weights==1)) {
	temp$weights <- weights
	temp$nevent.wt <- ff$event2[keep]
	temp$nrisk.wt  <- ff$nrisk2[keep]
	}
    temp
    }
