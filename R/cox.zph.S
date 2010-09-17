# $Id: cox.zph.S 11218 2009-02-09 12:09:29Z therneau $
#  Test proportional hazards
#
cox.zph <- function(fit, transform='km', global=TRUE) {
    call <- match.call()
    if (!inherits(fit, 'coxph')) stop ("Argument must be the result of coxph")
    if (inherits(fit, 'coxph.null'))
	stop("The are no score residuals for a Null model")

    sresid <- resid(fit, 'schoenfeld')
    varnames <- names(fit$coefficients)
    nvar <- length(varnames)
    ndead<- length(sresid)/nvar
    if (nvar==1) times <- as.numeric(names(sresid))
    else         times <- as.numeric(dimnames(sresid)[[1]])

    # Next line is no longer necessary: survfit.km can handle (start,stop] data
    # if (missing(transform) && attr(fit$y, 'type') != 'right')
    #      transform <- 'identity'
    if (is.character(transform)) {
	tname <- transform
	ttimes <- switch(transform,
			   'identity'= times,
			   'rank'    = rank(times),
			   'log'     = log(times),
			   'km' = {
				temp <- survfitKM(factor(rep(1,nrow(fit$y))),
						    fit$y, se.fit=FALSE)
				# A nuisance to do left cont KM
				t1 <- temp$surv[temp$n.event>0]
				t2 <- temp$n.event[temp$n.event>0]
				km <- rep(c(1,t1), c(t2,0))
				if (is.null(attr(sresid, 'strata')))
				    1-km
				else (1- km[sort.list(sort.list(times))])
				},
			   stop("Unrecognized transform"))
	}
    else {
	tname <- deparse(substitute(transform))
        if (length(tname) >1) tname <- 'user'
	ttimes <- transform(times)
	}
    xx <- ttimes - mean(ttimes)

    r2 <- sresid %*% fit$var * ndead
    test <- xx %*% r2        # time weighted col sums
    corel <- c(cor(xx, r2))
    z <- c(test^2 /(diag(fit$var)*ndead* sum(xx^2)))
    Z.ph <- cbind(corel, z, 1- pchisq(z,1))

    if (global && nvar>1) {
	test <- c(xx %*% sresid)
	z    <- c(test %*% fit$var %*% test) * ndead / sum(xx^2)
	Z.ph <- rbind(Z.ph, c(NA, z, 1-pchisq(z, ncol(sresid))))
	dimnames(Z.ph) <- list(c(varnames, "GLOBAL"), c("rho", "chisq", "p"))
	}
    else dimnames(Z.ph) <- list(varnames, c("rho", "chisq", "p"))

    dimnames(r2) <- list(times, names(fit$coefficients))
    temp <-list(table=Z.ph, x=ttimes, 
                y=r2 + outer(rep(1,ndead), fit$coefficients),
                var=fit$var, call=call, transform=tname)

    if (is.R()) class(temp) <- "cox.zph"
    else oldClass(temp) <- "cox.zph"
    temp
    }

"[.cox.zph" <- function(x, ..., drop=FALSE) {
    i <- ..1
    z<- list(table=x$table[i,,drop=FALSE], x=x$x, y=x$y[ ,i,drop=FALSE],
		var=x$var[i,i, drop=FALSE], call=x$call,
		transform=x$transform)
    attributes(z) <- attributes(x)
    z
    }
