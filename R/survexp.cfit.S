#
#
#  Do expected survival based on a Cox model
#   A fair bit of the setup work is identical to survfit.coxph, i.e.,
#     to reconstruct the data frame
#
#  The execution path for individual survival is completely separate, and
#    a whole lot simpler.
#
survexp.cfit <- function(group, x, y, death, individual, cox, se.fit, method) {
    if (!is.matrix(x)) stop("x must be a matrix")

    #
    # If it is individual survival, things are fairly easy
    #    (the parent routine has guarranteed NO strata in the Cox model
    #
    if (individual) {
	fit <- survfit.coxph(cox, se.fit=FALSE)
	risk <- x %*% cox$coefficients  -  
		           sum(cox$coefficients *cox$means)
	nt <- length(fit$time)
	surv <- approx(-c(0,fit$time), c(1,fit$surv), -y,
				method='constant', rule=2, f=1)$y
	return(list(times=y, surv=c(surv^(exp(risk)))))
	}

    # Otherwise, get on with the real work
    temp <- coxph.getdata(cox, y=TRUE, x=se.fit, stratax=FALSE)
    cy <- temp$y
    cx <- temp$x
    cn <- nrow(cy)
    nvar <- length(cox$coefficients)

    if (ncol(x) != nvar)
	stop("x matrix does not match the cox fit")

    ngrp <- group
    if (!is.logical(death)) stop("Invalid value for death indicator")

    if (missing(method))
	method <- (1 + 1*(cox$method=='breslow') +2*(cox$method=='efron')
		     + 10*(death))
    else stop("Program unfinished")

    #
    # Data appears ok so proceed
    #  First sort the old data set
    # Also, expand y to (start, stop] form.  This leads to slower processing,
    #  but I only have to program one case instead of 2.
    if (ncol(cy) ==2) {
  	mintime <- min(cy[,1])
 	if (mintime < 0) cy <- cbind(2*mintime-1, cy)
 	else	       cy <- cbind(-1, cy)
 	}
    ord <- order(cy[,2], -cy[,3])
    cy  <- cy[ord,]
    score <- exp(cox$linear.predictors[ord])
    if (se.fit) cx <- cx[ord,]
    else  cx <- 0   #dummy, for .C call


    #
    # Process the new data
    #
    if (missing(y) || is.null(y)) y <- rep(max(cy[,2]), nrow(x))
    ord <- order(group)
    group <- group - min(group)
    n <- nrow(x)
    ncurve <- length(unique(group))
    utimes <- unique(cy[cy[,3]==1,2])
    npt <- length(utimes)  #unique death times
    storage.mode(cy) <- 'double'
    xxx  <- .C(Cagsurv3, as.integer(n),
			  as.integer(nvar),
			  as.integer(ncurve),
			  as.integer(npt),
			  as.integer(se.fit),
			  as.double(score),
			  y = as.double(y[ord]),
                          as.integer(group[ord]),
			  x[ord,],
			  cox$coefficients,
			  cox$var,
			  cox$means,
			  as.integer(cn),
			  cy,
			  as.double(cx),
			  surv = matrix(0.0, npt, ncurve),
			  varhaz = matrix(0.0, npt, ncurve),
			  nrisk  = matrix(0.0, npt, ncurve),
			  as.integer(method), DUP=FALSE)

    surv <- apply(xxx$surv, 2, cumprod)
    if (se.fit)
	list(surv=surv, n=xxx$nrisk, times=utimes,
			se=sqrt(xxx$varhaz)/surv)
    else
	list(surv=surv, n=xxx$nrisk, times=utimes )
    }
