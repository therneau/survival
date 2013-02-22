# 
#  Do the actual fit of a survreg model.  This routine is for the case
#   of no penalized terms (splines, etc).
#
survreg.fit<- function(x, y, weights, offset, init, controlvals, dist, 
		       scale=0, nstrat=1, strata, parms=NULL, assign) {

    iter.max <- controlvals$iter.max
    eps <- controlvals$rel.tolerance
    toler.chol <- controlvals$toler.chol
    debug <- controlvals$debug

    if (!is.matrix(x)) stop("Invalid X matrix ")
    n <- nrow(x)
    nvar <- ncol(x)
    ny <- ncol(y)
    if (is.null(offset)) offset <- rep(0,n)
    if (missing(weights)|| is.null(weights)) weights<- rep(1.0,n)
    else if (any(weights<=0)) stop("Invalid weights, must be >0")

    if (scale <0) stop("Invalid scale")
    if (scale >0 && nstrat >1) 
	    stop("Cannot have both a fixed scale and strata")
    if (nstrat>1 && (missing(strata) || length(strata)!= n))
	    stop("Invalid strata variable")
    if (nstrat==1) strata <- rep(1,n)
    if (scale >0) nstrat2 <- 0          # number of variances to estimate
    else          nstrat2 <- nstrat

    if (is.character(dist)) {
	sd <- survreg.distributions[[dist]]
	if (is.null(sd)) stop ("Unrecognized distribution")
	}
    else sd <- dist
    if (!is.function(sd$density)) 
	stop("Missing density function in the definition of the distribution")
    dnum <- match(sd$name, c("Extreme value", "Logistic", "Gaussian"))
    if (is.na(dnum)) {
        #  We need to set up a callback routine
        #  This returns the 5 number distribution summary (see the density
        #  functions in survreg.distributions).  Interval censored obs require
        #  2 evals and all others 1, so the call to the routine will have n2
        #  values.
	dnum <- 4  # flag for the C routine
	n2 <- n + sum(y[,ny]==3)  
        
	#
        # Create an expression that will be evaluated by the C-code,
        #   but with knowledge of some current variables
        # In the R doc, this would be "body(function(z) {"
	#  in Splus (Chambers book):  "functionBody(function(z)"
	#  same action, different name.  Luckily 'quote' exists in both.
	# We make very sure the result is the right type and length here
	#  rather than in the C code, for simplicity.
        f.expr <- quote({
                if (length(parms)) temp <- sd$density(z, parms)
                else               temp <- sd$density(z)
                if (!is.matrix(temp) || any(dim(temp) != c(n2,5)) ||
                    !is.numeric(temp))
		    stop("Density function returned an invalid matrix")
                as.vector(as.double(temp))
                })
       
        # create an isolated sandbox (frame or environment) in which
	#  we can do the evaluation without endangering local objects
	#  but still with knowlege of sd, parms, and n2
        if (is.R()) rho <- new.env() #inherits necessary objects
	# SPlus else        rho <- new.frame(list(sd=sd, parms=parms, n2=n2))
        }
    else {
	f.expr <- 1  #dummy values for the .Call
	rho <- 1
	}

    # This is a subset of residuals.survreg: define the first and second
    #   derivatives at z=0 for the 4 censoring types
    #   Used below for starting estimates
    derfun <- function(y, eta, sigma, density, parms) {
	ny <- ncol(y)
	status <- y[,ny]
	z <- (y[,1] - eta)/sigma
	dmat <- density(z,parms)
	dtemp<- dmat[,3] * dmat[,4]    #f'
	if (any(status==3)) {
	    z2 <- (y[,2] - eta)/sigma
	    dmat2 <- density(z2, parms)
	    }
	else {
	    dmat2 <- matrix(0,1,5)   #dummy values
	    z2 <- 0
	    }
	tdenom <- ((status==0) * dmat[,2]) +
		  ((status==1) * 1 )       +
		  ((status==2) * dmat[,1]) +
		  ((status==3) * ifelse(z>0, dmat[,2]-dmat2[,2], 
		                             dmat2[,1] - dmat[,1]))
	tdenom <- 1/(tdenom* sigma)
	dg <- -tdenom   *(((status==0) * (0-dmat[,3])) +
			  ((status==1) * dmat[,4]) + 
			  ((status==2) * dmat[,3]) +
			  ((status==3) * (dmat2[,3]- dmat[,3])))

	ddg <- (tdenom/sigma)*(((status==0) * (0- dtemp)) +
			       ((status==1) * dmat[,5]) +
			       ((status==2) * dtemp) +
			       ((status==3) * (dmat2[,3]*dmat2[,4] - dtemp))) 
	list(dg = dg, ddg = ddg - dg^2)
	}

    #
    # A good initial value of the scale turns out to be critical for successful
    #   iteration, in a surprisingly large number of data sets.
    # The best way we've found to get one is to fit a model with only the
    #   mean and the scale.  We don't need to do this in 3 situations:
    #    1. The only covariate is a mean (this step is then just a duplicate
    #       of the main fit).
    #    2. There are no scale parameters to estimate
    #    3. The user gave initial estimates for the scale
    # However, for 2 and 3 we still want the loglik for a mean only model
    #  as a part of the returned object.
    #
    nvar2 <- nvar + nstrat2
    meanonly <- (nvar==1 && all(x==1))
    if (!meanonly) {
	yy <- ifelse(y[,ny]!=3, y[,1], (y[,1]+y[,2])/2 )
	coef <- sd$init(yy, weights, parms)  #starting estimate for this model
	#init returns \sigma^2, I need log(sigma)
	# We sometimes get into trouble with a small estimate of sigma,
	#  (the surface isn't SPD), but never with a large one.  Double it.
	if (scale >0) vars <- log(scale)
	else vars <- log(4*coef[2])/2    # log(2*sqrt(variance)) = log(4*var)/2
	coef <- c(coef[1], rep(vars, nstrat))
	
	# get a better initial value for the mean using the "glim" trick
	deriv <- derfun(y, yy, exp(vars), sd$density, parms)
	wt <-  -1*deriv$ddg*weights
	coef[1] <- sum(weights*deriv$dg + wt*(yy -offset)) / sum(wt)
        
	# Now the fit proper (intercept only)
	fit0 <- .Call(Csurvreg6,
		       iter = as.integer(20),
		       nvar = as.integer(1),
		       as.double(y),
		       as.integer(ny),
		       x = as.double(rep(1.0, n)),
		       as.double(weights),
		       as.double(offset),
		       coef= as.double(coef),
		       as.integer(nstrat2),
		       as.integer(strata),
		       as.double(eps),
		       as.double(toler.chol), 
		       as.integer(dnum),
		       f.expr,
		       rho)
	}
    #
    # Fit the model with all covariates
    #
    if (is.numeric(init)) {
	if (length(init) == nvar && (nvar2 > nvar)) {
	    # Add on the variance estimates from above
	    init <- c(init, fit0$coef[-1])
	    }
	if (length(init) != nvar2) stop("Wrong length for initial parameters")
	if (scale >0) init <- c(init, log(scale))
	}
    else  {
	# Do the 'glim' method of finding an initial value of coef
	if (meanonly) {
	    yy <- ifelse(y[,ny]!=3, y[,1], (y[,1]+y[,2])/2 )
	    coef <- sd$init(yy, weights, parms)
	    if (scale >0) vars <- rep(log(scale), nstrat)
	    else vars  <- rep(log(4*coef[2])/2, nstrat)  
	    }
	else vars <- fit0$coef[-1]
	eta <- yy - offset     #what would be true for a 'perfect' model

	deriv <- derfun(y, yy, exp(vars[strata]), sd$density, parms)
	wt <-  -1*deriv$ddg*weights
	coef <- coxph.wtest(t(x)%*% (wt*x), 
		       c((wt*eta + weights*deriv$dg)%*% x),
			    toler.chol=toler.chol)$solve
	init <- c(coef, vars)
	}

    # Now for the fit in earnest
    fit <- .Call(Csurvreg6,
		   iter = as.integer(iter.max),
		   as.integer(nvar),
		   as.double(y),
		   as.integer(ny),
		   as.double(x),
	           as.double(weights),
		   as.double(offset),
		   as.double(init),
	           as.integer(nstrat2),
	           as.integer(strata),
		   as.double(eps),
	           as.double(toler.chol), 
		   as.integer(dnum),
		   f.expr,
		   rho)

    if (iter.max >1 && fit$flag > nvar2) {
        warning("Ran out of iterations and did not converge")
	}

    cname <- dimnames(x)[[2]]
    if (is.null(cname)) cname <- paste("x", 1:ncol(x))
    if (scale==0) cname <- c(cname, rep("Log(scale)", nstrat))
    if (scale>0) fit$coef <- fit$coef[1:nvar2]
    names(fit$coef) <- cname

    if (meanonly) {
	coef0 <- fit$coef
	loglik <- rep(fit$loglik,2)
	}
    else {
	coef0 <- fit0$coef
	names(coef0) <- c("Intercept", rep("Log(scale)", nstrat))
	loglik <- c(fit0$loglik, fit$loglik)
	}
    temp <- list(coefficients   = fit$coef,
		 icoef  = coef0, 
		 var    = matrix(fit$var, nvar2, dimnames=list(cname, cname)),
		 loglik = loglik, 
		 iter   = fit$iter,
		 linear.predictors = c(x %*% fit$coef[1:nvar] + offset),
                 df= length(fit$coef),
		 score = fit$u
		 )
    temp
    }
