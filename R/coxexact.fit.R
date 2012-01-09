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
	newstrat <- as.integer(c(1*(diff(strata)!=0), 1))
	}
    y <- y[sorted,]
    x <- x[sorted,]
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
    newx <- scale(x)
    
    cfit <- .Call("Coxexact", 
                  as.integer(maxiter),
                  as.double(y),  # interger data?  Just in case.
                  newx,
                  as.double(offset),
                  as.integer(newstrat),
                  as.double(init),
                  as.double(control$eps),
                  as.double(control$toler.chol)
              )

    var <- matrix(agfit$imat,nvar,nvar)
    coef <- agfit$coef
    if (agfit$flag < nvar) which.sing <- diag(var)==0
    else which.sing <- rep(FALSE,nvar)

    infs <- abs(agfit$u %*% var)
    if (control$iter.max >1) {
	if (agfit$flag == 1000)
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
    lp  <- x %*% coef + offset - sum(coef *agfit$means)
    score <- as.double(exp(lp[sorted]))
    agres <- .C("agmart",
		   as.integer(n),
		   as.integer(0),
		   sstart, sstop,
		   sstat,
		   score,
		   rep(1.0, n),
		   newstrat,
		   resid=double(n),
                   PACKAGE = 'survival')
    resid <- double(n)
    resid[sorted] <- agres$resid
    names(resid) <- rownames
    coef[which.sing] <- NA

    list(coefficients  = coef,
		var    = var,
		loglik = agfit$loglik,
		score  = agfit$sctest,
		iter   = agfit$iter,
		linear.predictors = lp,
		residuals = resid,
		means = agfit$means,
		method= 'coxph')
    }
