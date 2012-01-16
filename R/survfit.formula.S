survfit.formula <- function(formula, data, weights, subset, 
			    na.action, etype, id, ...) {
    Call <- match.call()
    Call[[1]] <- as.name('survfit')  #make nicer printout for the user
    # create a copy of the call that has only the arguments we want,
    #  and use it to call model.frame()
    mfnames <- c('formula', 'data', 'weights', 'subset','na.action',
                 'etype', 'id')  #legal args for model.frame
    temp <- Call[c(1, match(mfnames, names(Call), nomatch=0))]
    temp[[1]] <- as.name("model.frame")
    if (is.R()) m <- eval.parent(temp)
    else        m <- eval(temp, sys.parent())
    
    Terms <- terms(formula, 'strata')
    ord <- attr(Terms, 'order')
    if (length(ord) & any(ord !=1))
	    stop("Interaction terms are not valid for this function")

    n <- nrow(m)
    Y <- model.extract(m, 'response')
    if (!is.Surv(Y)) stop("Response must be a survival object")

    casewt <- model.extract(m, "weights")
    if (is.null(casewt)) casewt <- rep(1,n)

    if (!is.null(attr(Terms, 'offset'))) warning("Offset term ignored")

    ll <- attr(Terms, 'term.labels')
    if (length(ll) == 0) X <- factor(rep(1,n))  # ~1 on the right
    else X <- strata(m[ll])
    
    etype <- model.extract(m, 'etype')
    id    <- model.extract(m, 'id')

    if (!is.Surv(Y)) stop("y must be a Surv object")

    # At one point there were lines here to round the survival
    # times to a certain number of digits.  This approach worked
    # almost all the time, but only almost.  The better logic is
    # now in the individual compuation routines
    if (attr(Y, 'type') != 'right' && attr(Y, 'type') != 'counting'){
        if (!is.null(etype)) 
            stop(paste("Cumulative incidence computation doesn't support \"", 
                       etype, "\" survival data", sep=''))
        else temp <-  survfitTurnbull(X, Y, casewt, ...)
    }
    else {
        # Sometimes a user will ask for a competing risk curve, but there is
        #  only one event type (often a result of a subset= clause).  In that
        #  case the usual KM is the same and faster.
        deaths <- Y[,ncol(Y)] ==1
        if (!is.null(etype) && length(unique(etype[deaths]))>1)
            temp <- survfitCI(X, Y, weights=casewt, etype=etype, id=id,  ...)
        else temp <- survfitKM(X, Y, casewt, ...)
    }
    
    class(temp) <- 'survfit'
    if (!is.null(attr(m, 'na.action')))
	    temp$na.action <- attr(m, 'na.action')

    temp$call <- Call
    temp
    }









