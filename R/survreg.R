survreg <- function(formula, data, weights, subset, na.action,
	dist='weibull', init=NULL,  scale=0, control, parms=NULL, 
	model=FALSE, x=FALSE, y=TRUE, robust=FALSE, cluster, score=FALSE,  ...) {

    Call <- match.call()    # save a copy of the call

    # Move any cluster() term out of the formula, and make it an argument
    #  instead.  This makes everything easier.  But, I can only do that with
    #  a local copy, doing otherwise messes up future use of update() on
    #  the model object for a user stuck in "+ cluster()" mode.
    if (missing(formula)) stop("a formula argument is required")
    
    ss <- c("cluster", "offset")
    Terms <- if (missing(data)) terms(formula, specials=ss) else
                 terms(formula, specials=ss, data=data)
    tcl <- attr(Terms, 'specials')$cluster
    if (length(tcl) > 1) stop("a formula cannot have multiple cluster terms")

    if (length(tcl) > 0) { # there is one
        # subscripting of formulas is broken at least through R 3.5, if the
        #  formula contains an offset.  Adding offset to the "specials" above
        #  is just a sneaky way to find out if one is present, then call
        #  reformulate ourselves.  tt is a correct index into the row labels
        #  of the factors attribute, tt+1 to the variables attribute (which is
        #  a list, so you have to skip the "list" call).  The term.labels attr
        #  contains neither the response nor the offset, but does contain the
        #  interactions, which we need.  
        factors <- attr(Terms, 'factors')
        if (any(factors[tcl,] >1)) stop("cluster() cannot be in an interaction")
        if (attr(Terms, "response") ==0)
            stop("formula must have a Surv response")
        # reformulate with the response option puts ` ` around Surv, which messes
        #  up evaluation, hence the fancy dance to replace a piece rather
        #  than recreate
        temp <- attr(Terms, "term.labels")
        oo <- attr(Terms, 'specials')$offset
        if (!is.null(oo)) {
            # add the offset to the set of labels
            ooterm <- rownames(factors)[oo]
            if (oo < tcl) temp <- c(ooterm, temp)
            else temp <- c(temp, ooterm)
        }
        if (is.null(Call$cluster))
            Call$cluster <- attr(Terms, "variables")[[1+tcl]][[2]]
        else warning("cluster appears both in a formula and as an argument, formula term ignored")
        formula[[3]]      <- reformulate(temp[1-tcl])[[2]]

        Call$formula <- formula
    }

    indx <- match(c("formula", "data", "weights", "subset", "na.action",
                    "cluster"),
                  names(Call), nomatch=0) 
    if (indx[1] ==0) stop("A formula argument is required")
    temp <- Call[c(1,indx)]  # only keep the arguments we wanted
    temp[[1L]] <- quote(stats::model.frame)   # change the function called

    special <- c("strata")
    temp$formula <- if(missing(data)) terms(formula, special)
                    else              terms(formula, special, data=data)
    m <- eval(temp, parent.frame())

    Terms <- attr(m, 'terms')
    weights <- model.extract(m, 'weights')
    Y <- model.extract(m, "response")
    if (!inherits(Y, "Surv")) stop("Response must be a survival object")
    type <- attr(Y, "type")
    if (type== 'counting') 
        stop ("start-stop type Surv objects are not supported")
    if (type=="mright" || type=="mcounting") 
        stop("multi-state survival is not supported")
   
    cluster <- model.extract(m, "cluster")
    if (length(cluster)) {
        if (missing(robust)) robust <- TRUE
        cluster <- as.numeric(as.factor(cluster))
        }
    else if (robust) cluster <- 1:nrow(Y)
   
    strats <- attr(Terms, "specials")$strata
    dropx <- NULL
    if (length(strats)) {
        temp <- untangle.specials(Terms, 'strata', 1)
        dropx <- temp$terms
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
        attr(newTerms, 'intercept') <- attr(Terms, 'intercept')
        }
    else               newTerms <- Terms
    X <- model.matrix(newTerms, m)
    assign <- lapply(attrassign(X, newTerms)[-1], function(x) x-1)
    xlevels <- .getXlevels(newTerms, m)
    contr.save <- attr(X, 'contrasts')
    
    if (!all(is.finite(X)))
        stop("data contains an infinite predictor")
        
    n <- nrow(X)
    nvar <- ncol(X)

    offset<- model.offset(m) # R returns NULL if no offset, Splus a zero
    if (length(offset)==0 || all(offset==0)) offset <- rep(0.,n)
    
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
    Ysave <- Y  # for use in the y component
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

        
        assign <- attrassign(X, newTerms)
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
    fit$terms <- Terms
    fit$contrasts <- contr.save
    if (length(xlevels)) fit$xlevels <- xlevels
    fit$means <- apply(X,2, mean)
    if (!is.null(weights)) fit$weights <- weights
    fit$call <- Call
    fit$dist <- dist
    if (model) fit$model <- m
    if (x)     fit$x <- X
    if (y)     fit$y <- Ysave
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
        else fit$var <- crossprod(residuals.survreg(fit, 'dfbeta'))
        if (!model) fit$model <- NULL  # take it back out
        }

    # set singular coefficients to NA
    #  this is purposely not done until the residuals, etc. are computed
    singular <- (diag(fit$var)==0)[1:length(fit$coef)]
    if (any(singular)) fit$coef <- ifelse(singular, NA, fit$coef)

    na.action <- attr(m, "na.action")
    if (length(na.action)) fit$na.action <- na.action

    if (any(pterms)) class(fit) <- c('survreg.penal', 'survreg')
    else	         class(fit) <- 'survreg'
    fit
    }
