if (!is.R()) setOldClass(c('survreg.penal', 'survreg'))

survreg <- function(formula, data, weights, subset, na.action,
	dist='weibull', init=NULL,  scale=0, control, parms=NULL, 
	model=FALSE, x=FALSE, y=TRUE, robust=FALSE, score=FALSE,  ...) {

    Call <- match.call()   # save a copy of the call
    indx <- match(c("formula", "data", "weights", "subset", "na.action"),
                  names(Call), nomatch=0) 
    if (indx[1] ==0) stop("A formula argument is required")
    temp <- Call[c(1,indx)]  # only keep the arguments we wanted
    temp[[1]] <- as.name('model.frame')  # change the function called

    special <- c("strata", "cluster")
    temp$formula <- if(missing(data)) terms(formula, special)
                    else              terms(formula, special, data=data)
    if (is.R()) m <- eval(temp, parent.frame())
    else        m <- eval(temp, sys.parent())

    Terms <- attr(m, 'terms')
    weights <- model.extract(m, 'weights')
    Y <- model.extract(m, "response")
    if (!inherits(Y, "Surv")) stop("Response must be a survival object")

    strats <- attr(Terms, "specials")$strata
    cluster<- attr(Terms, "specials")$cluster
    dropx <- NULL
    if (length(cluster)) {
        if (missing(robust)) robust <- TRUE
        tempc <- untangle.specials(Terms, 'cluster', 1:10)
        ord <- attr(Terms, 'order')[tempc$terms]
        if (any(ord>1)) stop ("Cluster can not be used in an interaction")
        cluster <- strata(m[,tempc$vars], shortlabel=TRUE)  #allow multiples
        dropx <- tempc$terms
        }
    if (length(strats)) {
        temp <- untangle.specials(Terms, 'strata', 1)
        dropx <- c(dropx, temp$terms)
        if (length(temp$vars)==1) strata.keep <- m[[temp$vars]]
        else strata.keep <- strata(m[,temp$vars], shortlabel=TRUE)
        strata <- as.numeric(strata.keep)
	nstrata <- max(strata)
        }
    else {
	nstrata <- 1
	strata <- 0
	}

    if (length(dropx)) {
        newTerms <- Terms[-dropx]
        # R (version 2.7.1) adds intercept=T anytime you drop something
        if (is.R()) attr(newTerms, 'intercept') <- attr(Terms, 'intercept')
        }
    else               newTerms <- Terms
    X <- model.matrix(newTerms, m)
    if (is.R()) {
	 assign <- lapply(attrassign(X, newTerms)[-1], function(x) x-1)
         xlevels <- .getXlevels(newTerms, m)
         contr.save <- attr(X, 'contrasts')
         }
    else {
        assign <- lapply(attr(X, 'assign')[-1], function(x) x -1)
        xvars <- as.character(attr(newTerms, 'variables'))
        xvars <- xvars[-attr(newTerms, 'response')]
        if (length(xvars) >0) {
                xlevels <- lapply(m[xvars], levels)
                xlevels <- xlevels[!unlist(lapply(xlevels, is.null))]
                if(length(xlevels) == 0)
                        xlevels <- NULL
                }
        else xlevels <- NULL
        contr.save <- attr(X, 'contrasts')
        }

    n <- nrow(X)
    nvar <- ncol(X)

    offset<- model.offset(m) # R returns NULL if no offset, Splus a zero
    if (length(offset)==0 || all(offset==0)) offset <- rep(0.,n)

    type <- attr(Y, "type")
    if (type== 'counting') stop ("Invalid survival type")
    
    # The user can either give a distribution name, in which the distribution
    #   is found in the object survreg.distributions, or include a list object
    #   of the same format as is found there.
    if (is.character(dist)) {
        # partial matching of names in [[ is on its way out in R, so
        #   first use match.arg, e.g. turn 'exp' into 'exponential'
        dist <- match.arg(dist, names(survreg.distributions))
	dlist <- survreg.distributions[[dist]]
	if (is.null(dlist)) stop(paste(dist, ": distribution not found"))
	}
    else if (is.list(dist)) dlist <- dist
    else stop("Invalid distribution object")

    #
    #   Make sure it is legal
    #
    if (!survregDtest(dlist)) stop("Invalid distribution object")

    # If the distribution is a transformation of another, perform
    #   said transform.
    #  
    logcorrect <- 0   #correction to the loglik due to transformations
    if (!is.null(dlist$trans)) {
	tranfun <- dlist$trans
	exactsurv <- Y[,ncol(Y)] ==1
	if (any(exactsurv)) {
            if (is.null(weights))
                 logcorrect <- sum(log(dlist$dtrans(Y[exactsurv, 1])))
            else logcorrect <- sum(weights[exactsurv]*log(dlist$dtrans(Y[exactsurv, 1])))
         }
	if (type=='interval') {
	    if (any(Y[,3]==3))
		    Y <- cbind(tranfun(Y[,1:2]), Y[,3])
	    else Y <- cbind(tranfun(Y[,1]), Y[,3])
	    }
	else if (type=='left')
	     Y <- cbind(tranfun(Y[,1]), 2-Y[,2])
	else     Y <- cbind(tranfun(Y[,1]), Y[,2])
	if (!all(is.finite(Y))) 
	    stop("Invalid survival times for this distribution")
	}
    else {
	if (type=='left') Y[,2] <- 2- Y[,2]
	else if (type=='interval' && all(Y[,3]<3)) Y <- Y[,c(1,3)]
	}

    if (is.null(dlist$itrans)) itrans <- function(x) x
    else itrans <- dlist$itrans

    if (!is.null(dlist$scale)) {
        if (!missing(scale)) warning(paste(dlist$name, 
                           "has a fixed scale, user specified value ignored"))
        scale <- dlist$scale
        }

    if (!is.null(dlist$dist))
        if (is.atomic(dlist$dist)) dlist <- survreg.distributions[[dlist$dist]]
        else                       dlist <- dlist$dist
    
    # check for parameters
    ptemp <- dlist$parms
    if (is.null(ptemp)) {
        if (!is.null(parms)) stop(paste(dlist$name, 
                              "distribution has no optional parameters"))
        }
    else {
        if (!is.numeric(ptemp)) 
            stop("Default parameters must be a numeric vector")
        if (!missing(parms)) {
            temp <- unlist(parms)  # just in case they gave a list object
            indx <- match(names(temp), names(ptemp))
            if (any(is.na(indx))) stop("Invalid parameter names")
            ptemp[names(ptemp)] <- temp
            }
        parms <- ptemp
        }

    # An idea originally from Brian R: if the user gave a list of
    #  control values, use it, but if they did not give an explicit control
    #  argument assume that they mistakenly wrote control parameters as a
    #  part of the "..." or other arguments
    if (missing(control)) control <- survreg.control(...)
    else control <- do.call('survreg.control', control)

    # The any() construction below is to catch a user that mistakenly
    #  thinks that 'scale' can be used in a model with multiple strata, and
    #  so provided a vector of scale values.
    # (A 'perhaps should be be added someday' feature).
    if (any(scale < 0)) stop("Invalid scale value")
    if (any(scale >0) && nstrata >1) 
	    stop("The scale argument is not valid with multiple strata")

    # Check for penalized terms
    pterms <- sapply(m, inherits, 'coxph.penalty')
    if (any(pterms)) {
	pattr <- lapply(m[pterms], attributes)
	# 
	# the 'order' attribute has the same components as 'term.labels'
	#   pterms always has 1 more (response), sometimes 2 (offset)
	# drop the extra parts from pterms
	temp <- c(attr(Terms, 'response'), attr(Terms, 'offset'))
	if (length(dropx)) temp <- c(temp, dropx+1)
	pterms <- pterms[-temp]
	temp <- match((names(pterms))[pterms], attr(Terms, 'term.labels'))
	ord <- attr(Terms, 'order')[temp]
	if (any(ord>1)) stop ('Penalty terms cannot be in an interaction')

        
        if (is.R()) assign <- attrassign(X, newTerms)
        else        assign <- attr( X, 'assign')   
        pcols <- assign[match(names(pterms[pterms]), names(assign))] 

        fit <- survpenal.fit(X, Y, weights, offset, init=init,
				controlvals = control,
			        dist= dlist, scale=scale,
			        strata=strata, nstrat=nstrata,
				pcols, pattr, parms=parms, assign)
	}
    else fit <- survreg.fit(X, Y, weights, offset, 
			    init=init, controlvals=control,
			    dist= dlist, scale=scale, nstrat=nstrata, 
			    strata, parms=parms)

    if (is.character(fit))  fit <- list(fail=fit)  #error message
    else {
	if (scale==0) {
	    nvar <- length(fit$coefficients) - nstrata
	    fit$scale <- exp(fit$coefficients[-(1:nvar)])
	    if (nstrata==1) names(fit$scale) <- NULL
	    else names(fit$scale) <- levels(strata.keep)
	    fit$coefficients  <- fit$coefficients[1:nvar]
	    fit$idf  <- 1 + nstrata
	    }
	else {
	    fit$scale <- scale
	    fit$idf  <- 1
	    }
	fit$loglik <- fit$loglik + logcorrect
	}

    if (!score) fit$score <- NULL   #do not return the score vector
    fit$df.residual <- n - sum(fit$df)
#   fit$fitted.values <- itrans(fit$linear.predictors)
    fit$terms <- Terms
    fit$contrasts <- contr.save
    if (length(xlevels)) fit$xlevels <- xlevels
    fit$means <- apply(X,2, mean)
    fit$call <- Call
    fit$dist <- dist
    if (model) fit$model <- m
    if (x)     fit$x <- X
    if (y)     fit$y <- Y
    if (length(parms)) fit$parms <- parms

    # Do this before attaching the na.action, so that residuals() won't
    #   reinsert missing values under na.exclude
    if (robust) {
        fit$naive.var <- fit$var
        if (!model) fit$model <- m  #temporary addition, so resid doesn't
                                    # have to reconstruct
        if (length(cluster))
             fit$var <- crossprod(rowsum(residuals.survreg(fit, 'dfbeta'), 
                                         cluster))
        else fit$var <- crossprod(rowsum(residuals.survreg(fit, 'dfbeta')))
        if (!model) fit$model <- NULL  # take it back out
        }

    na.action <- attr(m, "na.action")
    if (length(na.action)) fit$na.action <- na.action

    if (is.R()) {
        if (any(pterms)) class(fit) <- c('survreg.penal', 'survreg')
        else	         class(fit) <- 'survreg'
        }
    else {
        if (any(pterms)) oldClass(fit) <- 'survreg.penal'
        else	     oldClass(fit) <- 'survreg'
        }
    fit
    }
