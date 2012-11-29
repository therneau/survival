survfit.formula <- function(formula, data, weights, subset, 
			    na.action, etype, id, istate, ...) {
    Call <- match.call()
    Call[[1]] <- as.name('survfit')  #make nicer printout for the user
    # create a copy of the call that has only the arguments we want,
    #  and use it to call model.frame()
    mfnames <- c('formula', 'data', 'weights', 'subset','na.action',
                 'istate', 'id', "etype")  #legal args for model.frame
    temp <- Call[c(1, match(mfnames, names(Call), nomatch=0))]
    temp[[1]] <- as.name("model.frame")
    if (is.R()) m <- eval.parent(temp)
    else        m <- eval(temp, sys.parent())
    
    Terms <- terms(formula, c("strata", "cluster"))
    ord <- attr(Terms, 'order')
    if (length(ord) & any(ord !=1))
	    stop("Interaction terms are not valid for this function")

    n <- nrow(m)
    Y <- model.extract(m, 'response')
    if (!is.Surv(Y)) stop("Response must be a survival object")

    casewt <- model.extract(m, "weights")
    if (is.null(casewt)) casewt <- rep(1,n)

    if (!is.null(attr(Terms, 'offset'))) warning("Offset term ignored")

    id    <- model.extract(m, 'id')
    istate <- model.extract(m,"istate")
    temp <- untangle.specials(Terms, "cluster")
    if (length(temp$vars)>0) {
        if (length(temp$vars) > 1) stop("can not have two cluster terms")
        if (!is.null(id)) stop("can not have both a cluster term and an id variable")       
        id <- m[[temp$vars]]
        Terms <- Terms[-temp$terms]
    }

    ll <- attr(Terms, 'term.labels')
    if (length(ll) == 0) X <- factor(rep(1,n))  # ~1 on the right
    else X <- strata(m[ll])
    
    if (!is.Surv(Y)) stop("y must be a Surv object")
    
    # Backwards support for the now-depreciated etype argument
    etype <- model.extract(m, "etype")
    if (!is.null(etype)) {
        if (attr(Y, "type") == "mcounting" ||
            attr(Y, "type") == "mright")
            stop("cannot use both the etype argument and mstate survival type")
        if (length(istate)) 
            stop("cannot use both the etype and istate arguments")
        status <- Y[,ncol(Y)]
        etype <- as.factor(etype)
        temp <- table(etype, status==0)

        if (all(rowSums(temp==0) ==1)) {
            # The user had a unique level of etype for the censors
            newlev <- levels(etype)[order(-temp[,2])] #censors first
        }
        else newlev <- c(" ", levels(etype)[temp[,1] >0])
        status <- factor(ifelse(status==0,0, as.numeric(etype)),
                             labels=newlev)

        if (attr(Y, 'type') == "right")
            Y <- Surv(Y[,1], status, type="mstate")
        else if (attr(Y, "type") == "counting")
            Y <- Surv(Y[,1], Y[,2], status, type="mstate")
        else stop("etype argument incompatable with survival type")
    }
                         
    # At one point there were lines here to round the survival
    # times to a certain number of digits.  This approach worked
    # almost all the time, but only almost.  The better logic is
    # now in the individual compuation routines
    if (attr(Y, 'type') == 'left' || attr(Y, 'type') == 'interval')
        temp <-  survfitTurnbull(X, Y, casewt, ...)
    else if (attr(Y, 'type') == "right" || attr(Y, 'type')== "counting")
        temp <- survfitKM(X, Y, casewt, ...)
    else if (attr(Y, 'type') == "mright" || attr(Y, "type")== "mcounting")
        temp <- survfitCI(X, Y, weights=casewt, id=id,  istate=istate, ...)
    else {
        # This should never happen
        stop("unrecognized survival type")
    }
    if (is.null(temp$states)) class(temp) <- 'survfit'
    else class(temp) <- c("survfitms", "survfit")

    if (!is.null(attr(m, 'na.action')))
	    temp$na.action <- attr(m, 'na.action')

    temp$call <- Call
    temp
    }









