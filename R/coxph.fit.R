coxph.fit <- function(x, y, strata, offset, init, control,
                      weights, method, rownames, resid=TRUE, 
                      nocenter=NULL)
{
    n <-  nrow(y)
    if (is.matrix(x)) nvar <- ncol(x)
    else {
	if (length(x)==0) nvar <-0
	else nvar <-1
	}
    time <- y[,1]
    status <- y[,2]

    # Sort the data (or rather, get a list of sorted indices)
    if (length(strata)==0) {
	sorted <- order(time)
        strata <- NULL
	newstrat <- as.integer(rep(0,n))
	}
    else {
	sorted <- order(strata, time)
	strata <- strata[sorted]
	newstrat <- as.integer(c(1*(diff(as.numeric(strata))!=0), 1))
	}
    if (missing(offset) || is.null(offset)) offset <- rep(0,n)
    if (missing(weights)|| is.null(weights))weights<- rep(1,n)
    else {
	if (any(weights<=0)) stop("Invalid weights, must be >0")
	weights <- weights[sorted]
	}
    stime <- as.double(time[sorted])
    sstat <- as.integer(status[sorted])

    if (nvar==0) {
	# A special case: Null model.
	#  (This is why I need the rownames arg- can't use x' names)
	# Set things up for 0 iterations on a dummy variable
	x <- as.matrix(rep(1.0, n))
	nullmodel <- TRUE
	nvar <- 1
	init <- 0
	maxiter <- 0
	}
    else {
	nullmodel <- FALSE
	maxiter <- control$iter.max
	if (!missing(init) && length(init)>0) {
	    if (length(init) != nvar) stop("Wrong length for inital values")
	    }
	else init <- rep(0,nvar)
	}

    # 2012 change: individually choose which variable to rescale
    # default: leave 0/1 variables alone
    if (is.null(nocenter)) zero.one <- rep(FALSE, ncol(x))
    else zero.one <- apply(x, 2, function(z) all(z %in% nocenter))

    storage.mode(weights) <- storage.mode(init) <- "double"
    coxfit <- .Call(Ccoxfit6, 
                     as.integer(maxiter),
                     stime, 
                     sstat,
                     x[sorted,],
                     as.double(offset[sorted]),
                     weights,
                     newstrat,
                     as.integer(method=="efron"),
                     as.double(control$eps),
                     as.double(control$toler.chol),
                     as.vector(init),
                     ifelse(zero.one, 0L, 1L))

    if (nullmodel) {
        if (resid) {
            score <- exp(offset[sorted])
            coxres <- .C(Ccoxmart, as.integer(n),
				as.integer(method=='efron'),
				stime,
				sstat,
				newstrat,
				as.double(score),
				as.double(weights),
				resid=double(n))
            resid <- double(n)
            resid[sorted] <- coxres$resid
            names(resid) <- rownames

            list( loglik = coxfit$loglik[1],
                 linear.predictors = offset,
                 residuals = resid,
                 method = method,
                 class = c('coxph.null', 'coxph') )
        }
	else list(loglik = coxfit$loglik[1],
                 linear.predictors = offset,
                 method = method,
                 class = c('coxph.null', 'coxph') )
    }
    else {
	coef <- coxfit$coef
	lp <- c(x %*% coef) + offset - sum(coef*coxfit$means)
	var <- matrix(coxfit$imat,nvar,nvar)
	if (coxfit$flag < nvar) which.sing <- diag(var)==0
	else which.sing <- rep(FALSE,nvar)

	infs <- abs(coxfit$u %*% var)
	if (maxiter >1) {
	    if (coxfit$flag == 1000) {
		   warning("Ran out of iterations and did not converge")
                   if (max(lp) > 500 || any(!is.finite(infs)))
                       warning("one or more coefficients may be infinite")
               }
	    else {
		infs <- (!is.finite(coxfit$u) |
                         ((infs > control$eps) & 
			 infs > control$toler.inf*abs(coef)))
		if (any(infs))
		warning(gettextf("Loglik converged before variable %s; coefficient may be infinite.",
			  paste((1:nvar)[infs],collapse=",")))
		}
	    }

	names(coef) <- dimnames(x)[[2]]
        if (maxiter > 0) coef[which.sing] <- NA  #leave it be if iter=0 is set
        rval <- list(coefficients  = coef,
		    var    = var,
		    loglik = coxfit$loglik,
		    score  = coxfit$sctest,
		    iter   = coxfit$iter,
		    linear.predictors = as.vector(lp),
		    means = coxfit$means,
                    method = method,
 		    class ='coxph')
	if (resid) {
            temp <- lp[sorted]
            if (any(temp > log(.Machine$double.xmax))) {
                # prevent a failure message due to overflow
                #  this occurs with near-infinite coefficients
                temp <- temp + log(.Machine$double.xmax) - (1 + max(temp))
            }
            score <- exp(temp)
            coxres <- .C(Ccoxmart, as.integer(n),
				as.integer(method=='efron'),
				stime,
				sstat,
				newstrat,
				as.double(score),
				as.double(weights),
				resid=double(n))
            resid <- double(n)
            resid[sorted] <- coxres$resid
            names(resid) <- rownames
            # rval$residuals <- resid  # not good enough
            # some other code expects the ORDER of components in a coxph
            # object will never change: residuals right before means
            rval <- c(rval[1:6], list(residuals=resid), rval[-(1:6)])
        }
        rval
    }
}          

