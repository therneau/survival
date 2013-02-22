# This routine fits right censored data when the method is 
#  "exact".  The most common use for this option is matched
#   case-control data.
coxexact.fit <- function(x, y, strata, offset, init, control,
			  weights, method, rownames)
    {
    if (!is.matrix(x)) stop("Invalid formula for cox fitting function")
    if (!is.null(weights) && any(weights!=1))
	  stop("Case weights are not supported for the exact method")
    n <- nrow(x)
    nvar <- ncol(x)

    # The risk set addition in the C-code, which is the critically slow
    # part of the calculations, expects to have the data in sorted order:
    #   (large to small times) within strata
    if (length(strata)==0) {
	sorted <- order(-y[,1])
	newstrat <- as.integer(rep(0,n))
	}
    else {
	sorted <- order(strata, -y[,1])
	strata <- (as.numeric(strata))[sorted]
	newstrat <- as.integer(c(1, 1*(diff(strata)!=0)))
	}
    y <- y[sorted,]
    if (is.null(offset)) offset <- rep(0.,n)
    else offset <- offset[sorted]

    if (is.null(nvar)) {
	# A special case: Null model.  Not worth coding up
	stop("Cannot handle a null model + exact calculation (yet)")
	}

    if (!is.null(init)) {
	if (length(init) != nvar) stop("Wrong length for inital values")
	}
    else init <- rep(0.,nvar)

    # Prescale the data set to improve numerical accuracy.
    #  We will undo the scaling before finishing up.
    newx <- scale(x[sorted,])
#    newx <- scale(x, scale=NULL)  #debug
    rescale <- attr(newx, "scaled:scale")
    means   <- attr(newx, "scaled:center")

    cfit <- .Call(Ccoxexact, 
                  as.integer(control$iter.max),
                  as.double(y),  # interger data?  Just in case.
                  newx,
                  as.double(offset),
                  as.integer(newstrat),
                  as.double(init*rescale),
                  as.double(control$eps),
                  as.double(control$toler.chol)
              )

    loglik <- cfit$loglik[1:2]  #these are packed into one vector
    sctest <- cfit$loglik[3]
    iter <- cfit$loglik[5]
    flag <- cfit$loglik[4]
    var <- matrix(cfit$imat,nvar,nvar)
    coef <- cfit$coef

    if (flag < nvar) which.sing <- diag(var)==0
    else which.sing <- rep(FALSE,nvar)

    infs <- abs(cfit$u %*% var)
    if (control$iter.max >1) {
	if (flag == 1000)
	       warning("Ran out of iterations and did not converge")
	    else {
		infs <- ((infs > control$eps) & 
			 infs > control$toler.inf*abs(coef))
		if (any(infs))
		warning(paste("Loglik converged before variable ",
			  paste((1:nvar)[infs],collapse=","),
			  "; beta may be infinite. "))
		}
	}

    names(coef) <- dimnames(x)[[2]]
    lp  <- newx %*% coef + offset 
    score <- as.double(exp(lp))

    # Compute the residuals
    cxres <- .C(Ccoxmart2,
		   as.integer(n),
		   as.double(y[,1]),
		   as.integer(y[,2]),
                   newstrat,
                   score,
                   rep(1.0, n),  #weights
                   resid=double(n),
                   DUP=FALSE)
    resid <- double(n)
    resid[sorted] <- cxres$resid
    names(resid) <- rownames
    coef[which.sing] <- NA
    lp.unsort <- double(n)
    lp.unsort[sorted] <- lp

    scmat <- diag(1/rescale, nvar,nvar)
    list(coefficients  = coef/rescale,
		var    = scmat %*% var %*% scmat,
		loglik = loglik,
		score  = sctest,
		iter   = iter,
		linear.predictors = lp.unsort,
		residuals = resid,
		means = means,
		method= 'coxph')
    }
