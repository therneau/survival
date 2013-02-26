residuals.coxph <-
  function(object, type=c("martingale", "deviance", "score", "schoenfeld",
			  "dfbeta", "dfbetas", "scaledsch","partial"),
	    collapse=FALSE, weighted=FALSE, ...)
    {
    type <- match.arg(type)
    otype <- type
    if (type=='dfbeta' || type=='dfbetas') {
	type <- 'score'
	if (missing(weighted))
            weighted <- TRUE  # different default for this case
	}
    if (type=='scaledsch') type<-'schoenfeld'

    n <- length(object$residuals)
    rr <- object$residuals
    y <- object$y
    x <- object[['x']]  # avoid matching object$xlevels
    vv <- object$naive.var
    if (is.null(vv)) vv <- object$var
    weights <- object$weights
    if (is.null(weights)) weights <- rep(1,n)
    strat <- object$strata
    method <- object$method
    if (method=='exact' && (type=='score' || type=='schoenfeld'))
	stop(paste(type, 'residuals are not available for the exact method'))

    if (type == 'martingale' || type == 'partial')
        rr <- object$residuals
    else {
	# I need Y, and perhaps the X matrix (and strata)
	Terms <- object$terms
	if (!inherits(Terms, 'terms'))
		stop("invalid terms component of object")
	strats <- attr(Terms, "specials")$strata
	if (is.null(y)  ||  (is.null(x) && type!= 'deviance')) {
	    temp <- coxph.getdata(object, y=TRUE, x=TRUE, stratax=TRUE)
	    y <- temp$y
	    x <- temp$x
	    if (length(strats)) strat <- temp$strata
	    }

	ny <- ncol(y)
	status <- y[,ny,drop=TRUE]

	if (type != 'deviance') {
	    nstrat <- as.numeric(strat)
	    nvar <- ncol(x)
	    if (is.null(strat)) {
		ord <- order(y[,ny-1], -status)
		newstrat <- rep(0,n)
		}
	    else {
		ord <- order(nstrat, y[,ny-1], -status)
		newstrat <- c(diff(as.numeric(nstrat[ord]))!=0 ,1)
		}
	    newstrat[n] <- 1

	    # sort the data
	    x <- x[ord,]
	    y <- y[ord,]
	    score <- exp(object$linear.predictors)[ord]
	    }
	}

    #
    # Now I have gotton the data that I need-- do the work
    #
    if (type=='schoenfeld') {
  	if (ny==2) {
 	    mintime <- min(y[,1])
 	    if (mintime < 0) y <- cbind(2*mintime -1, y)
 	    else             y <- cbind(-1,y)
 	    }
	temp <- .C(Ccoxscho, n=as.integer(n),
			    as.integer(nvar),
			    as.double(y),
			    resid=  as.double(x),
			    as.double(score * weights[ord]),
			    as.integer(newstrat),
			    as.integer(method=='efron'),
			    double(3*nvar) )

	deaths <- y[,3]==1

	if (nvar==1) rr <- temp$resid[deaths]
	else         rr <- matrix(temp$resid[deaths], ncol=nvar) #pick rows 
	if (weighted) rr <- rr * weights[deaths]

	if (length(strats)) attr(rr, "strata")  <- table((strat[ord])[deaths])
	time <- c(y[deaths,2])  # 'c' kills all of the attributes
	if (is.matrix(rr)) dimnames(rr)<- list(time, names(object$coefficients))
	else               names(rr) <- time

	if (otype=='scaledsch') {
	    ndead <- sum(deaths)
	    coef <- ifelse(is.na(object$coefficients), 0, object$coefficients)
	    if (nvar==1) rr <- rr*vv *ndead + coef
	    else         rr <- rr %*%vv * ndead + outer(rep(1,nrow(rr)),coef)
	    }
	return(rr)
	}

    if (type=='score') {
	if (ny==2) {
	    resid <- .C(Ccoxscore, as.integer(n),
				as.integer(nvar),
				as.double(y),
				x=as.double(x),
				as.integer(newstrat),
				as.double(score),
				as.double(weights[ord]),
				as.integer(method=='efron'),
				resid= double(n*nvar),
				double(2*nvar))$resid
	    }
	else {
	    resid<- .C(Cagscore,
				as.integer(n),
				as.integer(nvar),
				as.double(y),
				as.double(x),
				as.integer(newstrat),
				as.double(score),
				as.double(weights[ord]),
				as.integer(method=='efron'),
				resid=double(n*nvar),
				double(nvar*6))$resid
	    }
	if (nvar >1) {
	    rr <- matrix(0, n, nvar)
	    rr[ord,] <- matrix(resid, ncol=nvar)
	    dimnames(rr) <- list(names(object$residuals), 
				 names(object$coefficients))
	    }
	else rr[ord] <- resid

	if      (otype=='dfbeta') {
	    if (is.matrix(rr)) rr <- rr %*% vv
	    else               rr <- rr * vv
	    }
	else if (otype=='dfbetas') {
	    if (is.matrix(rr))  rr <- (rr %*% vv) %*% diag(sqrt(1/diag(vv)))
	    else                rr <- rr * sqrt(vv)
	    }
	}
    
    #
    # Multiply up by case weights (which will be 1 for unweighted)
    #
    if (weighted) rr <- rr * weights
    
    #Expand out the missing values in the result
    if (!is.null(object$na.action)) {
	rr <- naresid(object$na.action, rr)
   	if (is.matrix(rr)) n <- nrow(rr)
	else               n <- length(rr)
	if (type=='deviance') status <- naresid(object$na.action, status)
	}
    
    if (type=="partial"){
	# This needs to be done after the naresid expansion, since the
	#   predict function will have done naresid expansion, so that
	#   the lengths match
        rr <- rr + predict(object,type="terms")
        }

    # Collapse if desired
    if (!missing(collapse)) {
	if (length(collapse) !=n) stop("Wrong length for 'collapse'")
	rr <- drop(rowsum(rr, collapse))
  	if (type=='deviance') status <- drop(rowsum(status, collapse))
	}

    # Deviance residuals are computed after collapsing occurs
    if (type=='deviance')
	sign(rr) *sqrt(-2* (rr+
			      ifelse(status==0, 0, status*log(status-rr))))
    else rr
    }
    
