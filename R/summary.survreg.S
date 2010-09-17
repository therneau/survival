# $Id: summary.survreg.S 11166 2008-11-24 22:10:34Z therneau $
summary.survreg<- function(object, correlation = FALSE,...)
{
    if (!is.null(object$fail)) {
	warning(" Survreg failed.", x$fail, "   No summary provided\n")
	return(invisible(object))
	}
    wt <- object$weights
    
    nvar0 <- length(object$coefficients)
    nvar <- nrow(object$var)
    if (nvar > nvar0) {
	coef <- c(object$coefficients, log(object$scale))
	if ( (nvar-nvar0)==1) cname <- c(names(object$coefficients), "Log(scale)")
	else cname <- c(names(object$coefficients), names(object$scale))
	}
    else {
	coef <- object$coefficients
	cname <- names(object$coefficients)
	}

    n <- length(object$linear.predictors)
    p <- sum(!is.na(coef))
    if(!p) {
        warning("This model has zero rank --- no summary is provided")
        return(invisible(object))
        }

    nsingular <- nvar - p
    if (is.null(object$naive.var)){
        table <- matrix(rep(coef, 4), ncol = 4)
        dimnames(table) <- list(cname, c("Value", "Std. Error", "z", "p"))
        stds <- sqrt(diag(object$var))
        table[, 2] <- stds
        table[, 3] <- table[, 1]/stds
        table[, 4] <- 2*pnorm(-abs(table[,3]))
    } else {
        table <- matrix(rep(coef, 5), ncol = 5)
        dimnames(table) <- list(cname, c("Value", "Std. Err","(Naive SE)", "z", "p"))
        stds <- sqrt(diag(object$var))
        table[, 2] <- stds
        table[, 3] <- sqrt(diag(object$naive.var))
        table[, 4] <- table[, 1]/stds
        table[, 5] <- 2*pnorm(-abs(table[,4]))
    }
    if(correlation) {
	nas <- is.na(coef)
	stds <- stds[!nas]
	correl <- diag(1/stds) %*% object$var[!nas, !nas] %*% diag(1/stds)
        dimnames(correl) <- list(cname, cname)
        }
    else correl <- NULL

    dist <- object$dist
    if (is.character(dist)) sd <- survreg.distributions[[dist]]
    else sd <- dist

    if (length(object$parms)) 
	    pprint<- paste(sd$name, 'distribution: parmameters=', object$parms)
    else    pprint<- paste(sd$name, 'distribution')

    x <- object[match(c('call', 'df', 'loglik', 'iter', 'na.action', 'idf',
			'scale', 'coefficients', 'var'), 
		      names(object), nomatch=0)]
    x <- c(x, list(table=table, correlation=correl, parms=pprint,
		   n=n, chi=2*diff(object$loglik)), 
	           robust=!is.null(object$naive.var))

    if (is.R()) class(x) <- 'summary.survreg'
    else        oldClass(x) <- 'summary.survreg'
    x
    }
