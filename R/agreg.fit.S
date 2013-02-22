agreg.fit <- function(x, y, strata, offset, init, control,
			weights, method, rownames)
    {
    n <- nrow(y)
    nvar <- ncol(x)
    start <- y[,1]
    stopp <- y[,2]
    event <- y[,3]
    if (all(event==0)) stop("Can't fit a Cox model with 0 failures")

    # Sort the data (or rather, get a list of sorted indices)
    #  For both stop and start times, the indices go from last to first
    if (length(strata)==0) {
	sort.end  <- order(-stopp, event)
	sort.start<- order(-start)
	newstrat  <- n
	}
    else {
	sort.end  <- order(strata, -stopp, event)
	sort.start<- order(strata, -start)
	newstrat  <- cumsum(table(strata))
	}
    if (missing(offset) || is.null(offset)) offset <- rep(0.0, n)
    if (missing(weights)|| is.null(weights))weights<- rep(1.0, n)
    else if (any(weights<=0)) stop("Invalid weights, must be >0")

    if (is.null(nvar) || nvar==0) {
	# A special case: Null model.  Just return obvious stuff
        #  To keep the C code to a small set, we call the usual routines, but
	#  with a dummy X matrix and 0 iterations
	nvar <- 1
	x <- matrix(1:n, ncol=1)  #keep the .C call happy
	maxiter <- 0
	nullmodel <- TRUE
        if (length(init) !=0) stop("Wrong length for inital values")
        init <- 0.0  #dummy value to keep the .C call happy (doesn't like 0 length)
	}
    else {
	nullmodel <- FALSE
	maxiter <- control$iter.max
        
        if (is.null(init)) init <- rep(0., nvar)
	if (length(init) != nvar) stop("Wrong length for inital values")
	}

    agfit <- .C(Cagfit3, iter= as.integer(maxiter),
		as.integer(n),
		as.integer(nvar), 
		as.double(start), 
		as.double(stopp),
		as.integer(event),
		as.double(x),
		as.double(offset),
		as.double(weights),
		as.integer(length(newstrat)),
		as.integer(newstrat),
		as.integer(sort.end-1),
		as.integer(sort.start-1),
		means = double(nvar),
		coef= as.double(init),
		u = double(nvar),
		imat= double(nvar*nvar), loglik=double(2),
		flag=integer(1),
		double(2*nvar*nvar +nvar*3 + n),
		as.double(control$eps),
		as.double(control$toler.chol),
		sctest=as.double(method=='efron'))

    var <- matrix(agfit$imat,nvar,nvar)
    coef <- agfit$coef
    if (agfit$flag < nvar) which.sing <- diag(var)==0
	else which.sing <- rep(FALSE,nvar)

    infs <- abs(agfit$u %*% var)
    if (maxiter >1) {
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
    lp  <- as.vector(x %*% coef + offset - sum(coef *agfit$means))
    score <- as.double(exp(lp))

    agres <- .C(Cagmart2,
		as.integer(n),
		as.integer(method=='efron'),
		as.double(start), 
		as.double(stopp),
		as.integer(event),
		as.integer(length(newstrat)), 
		as.integer(newstrat),
		as.integer(sort.end-1), 
		as.integer(sort.start-1),
		score,
		as.double(weights),
		resid=double(n),
		double(2*sum(event)))
    resid <- agres$resid

    if (nullmodel) {
	resid <- agres$resid
	names(resid) <- rownames

	list(loglik=agfit$loglik[2],
	     linear.predictors = offset,
	     residuals = resid,
	     method= c("coxph.null", 'coxph') )
	}
    else {
	names(coef) <- dimnames(x)[[2]]
	names(resid) <- rownames
	coef[which.sing] <- NA

        concordance <- survConcordance.fit(y, lp, strata, weights) 
	list(coefficients  = coef,
	     var    = var,
	     loglik = agfit$loglik,
	     score  = agfit$sctest,
	     iter   = agfit$iter,
	     linear.predictors = as.vector(lp),
	     residuals = resid,
	     means = agfit$means,
             concordance = concordance,
	     method= 'coxph')
	}
    }


