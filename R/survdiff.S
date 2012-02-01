survdiff <- function(formula, data, subset, na.action, rho=0) {
    call <- match.call()
    m <- match.call(expand.dots=FALSE)
    m$rho <- NULL

    if (!inherits(formula, 'formula'))
        stop("The 'formula' argument is not a formula")

    Terms <- if(missing(data)) terms(formula, 'strata')
	     else              terms(formula, 'strata', data=data)
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    if (is.R()) m <- eval(m, parent.frame())
    else        m <- eval(m, sys.parent())

    y <- model.extract(m, "response")
    if (!inherits(y, "Surv")) stop("Response must be a survival object")
    if (attr(y, 'type') != 'right') stop("Right censored data only")
    ny <- ncol(y)
    n <- nrow(y)

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
	    num <-  var - observed
	    }
	chi <- num*num/var
	rval <-list(n= n, obs = observed, exp=expected, var=var,
			chisq= chi)
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

	rval <-list(n= table(groups), obs = fit$observed,
		    exp = fit$expected, var=fit$var,  chisq=chi)
	if (length(strats)) rval$strata <- table(strata.keep)
	}

    na.action <- attr(m, "na.action")
    if (length(na.action)) rval$na.action <- na.action
    rval$call <- call

    if (is.R()) class(rval) <- 'survdiff'
    else        oldClass(rval) <- 'survdiff'
    rval
    }
