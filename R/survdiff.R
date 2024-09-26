survdiff <- function(formula, data, subset, na.action, rho=0, timefix=TRUE) {
    call <- match.call()
    m <- match.call(expand.dots=FALSE)
    m$rho <- NULL

    if (!inherits(formula, 'formula'))
        stop("The 'formula' argument is not a formula")
    formula <- removeDoubleColonSurv(formula)  # protect from survival::strata
    call$formula <- formula

    Terms <- if(missing(data)) terms(formula, 'strata')
	     else              terms(formula, 'strata', data=data)
    m$formula <- Terms
    # Make "strata" be local to the formula, without invoking any
    #  outside functions. We do this by inserting another environment on
    #  the front of the search path.  This is
    #  part of my defense against use of survival::strata.  Putting a local
    #  copy first on the path allows for users who don't want to load the
    #  survival namespace.
    coxenv <- new.env(parent= environment(formula))
    assign("strata", survival::strata, envir= coxenv)
    environment(m$formula) <- coxenv

    m[[1L]] <- quote(stats::model.frame)
    m <- eval(m, parent.frame())
    
    y <- model.extract(m, "response")
    if (!inherits(y, "Surv")) stop("Response must be a survival object")
    if (attr(y, "type") %in% c("mright", "mcounting"))
        stop("survdiff not defined for multi-state data")
    if (attr(y, "type") == "counting")
        stop("survdiff not defined for counting process data")
    if (attr(y, 'type') != "right") 
        stop("Right censored data only")
    ny <- ncol(y)
    n <- nrow(y)

     # Deal with the near-ties problem
    if (!is.logical(timefix) || length(timefix) > 1)
        stop("invalid value for timefix option")
    if (timefix) y <- aeqSurv(y) 
 
    offset<- attr(Terms, "offset")
    if (!is.null(offset)) {
	#one sample test
	offset <- as.numeric(m[[offset]])
	if (length(attr(Terms,"factors"))>0) 
		stop("Cannot have both an offset and groups")
	if (any(offset <0 | offset >1))
	    stop("The offset must be a survival probability")

	expected <- sum(-log(offset))  #sum of expected events
	observed <- sum(y[,ny])
	if (rho!=0) {
	    num <- sum(1/rho - ((1/rho + y[,ny])*offset^rho))
	    var <- sum(1- offset^(2*rho))/(2*rho)
	    }
	else {
	    var <-  sum(-log(offset))
	    num <-  observed - var
	    }
	chi <- num*num/var
	rval <-list(n= n, obs = observed, exp=expected, var=var,
		chisq= chi, pvalue = pchisq(chi, df=1, lower.tail=FALSE))
	}

    else { #k sample test
	strats <- attr(Terms, "specials")$strata
	if (length(strats)) {
	    temp <- untangle.specials(Terms, 'strata', 1)
	    dropx <- temp$terms
	    if (length(temp$vars)==1) strata.keep <- m[[temp$vars]]
	    else strata.keep <- strata(m[,temp$vars], shortlabel=TRUE)
	    }
	else strata.keep <- rep(1,nrow(m))

	#Now create the group variable
	if (length(strats)) ll <- attr(Terms[-dropx], 'term.labels')
	else                ll <- attr(Terms, 'term.labels')
	if (length(ll) == 0) stop("No groups to test")
	else groups <- strata(m[ll])

	fit <- survdiff.fit(y, groups, strata.keep, rho)
	if (is.matrix(fit$observed)){
	    otmp <- apply(fit$observed,1,sum)
	    etmp <- apply(fit$expected,1,sum)
	    }
	else {
	    otmp <- fit$observed
	    etmp <- fit$expected
	    }
	df   <- (etmp >0)                #remove groups with exp=0
	if (sum(df) <2) chi <- 0         # No test, actually
	else {
	    temp2 <- ((otmp - etmp)[df])[-1]
	    vv <- (fit$var[df,df])[-1,-1, drop=FALSE]
	    chi <- sum(solve(vv, temp2) * temp2)
	    }
	df <- (sum(1*(etmp>0))) -1
	rval <-list(n= table(groups), obs = fit$observed,
		    exp = fit$expected, var=fit$var,  chisq=chi, 
                    pvalue= pchisq(chi, df, lower.tail=FALSE))
	if (length(strats)) rval$strata <- table(strata.keep)
	}

    na.action <- attr(m, "na.action")
    if (length(na.action)) rval$na.action <- na.action
    rval$call <- call

    class(rval) <- 'survdiff'
    rval
    }
