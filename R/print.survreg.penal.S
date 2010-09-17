# $Id: print.survreg.penal.S 11166 2008-11-24 22:10:34Z therneau $
print.survreg.penal <-  function(x, terms=FALSE, maxlabel=25,
			       digits=max(options()$digits - 4, 3), ...) {
    if (!inherits(x, 'survreg.penal')) stop("Invalid object")

    if (!is.null(x$call)) {
	cat("Call:\n")
	dput(x$call) 
	cat("\n")
	}
    if (!is.null(x$fail)) {
	cat(" Survreg failed.", x$fail, "\n")
	return()
	}
    savedig <- options(digits = digits)
    on.exit(options(savedig))

    coef <- x$coefficients
    if (length(coef)==0)
	    stop("Penalized fits must have an intercept!")

    #
    # Map terms to special print functions, and the list of iteration histories
    #
    pterms <- x$pterms
    nterms <- length(pterms)
    npenal <- sum(pterms>0)
    print.map <- rep(0,nterms)
    if (!is.null(x$printfun)) {
	temp <- unlist(lapply(x$printfun, is.null))  #which ones are missing
	print.map[pterms>0] <- (1:npenal) * (!temp)
	}

    # Tedious, but build up the coef matrix a term at a time
    print1 <- NULL
    pname1 <- NULL
    if (is.null(x$assign2)) alist <- x$assign
    else alist <- x$assign2

    print2 <- NULL
    for (i in 1:nterms) {
	kk <- alist[[i]]
	if (print.map[i] >0) {
	    j <- print.map[i]	
	    if (pterms[i]==2) 
		 temp <- (x$printfun[[j]])(x$frail, x$fvar, ,x$df[i], 
					   x$history[[j]])
	    else temp <- (x$printfun[[j]])(coef[kk], x$var[kk,kk], 
					   x$var2[kk,kk], 
					   x$df[i], x$history[[j]])
	    print1 <- rbind(print1, temp$coef)
	    if (is.matrix(temp$coef)) {
		xx <- dimnames(temp$coef)[[1]]
		if (is.null(xx))
			xx <- rep(names(pterms)[i], nrow(temp$coef))
		else    xx <- paste(names(pterms)[i], xx, sep=', ')
		pname1 <- c(pname1, xx)
		}
	    else  pname1 <- c(pname1, names(pterms)[i])
	    print2 <- c(print2, temp$history)
	    }

	else if (terms && length(kk)>1) {
	    pname1 <- c(pname1, names(pterms)[i])
	    temp <- coxph.wtest(x$var[kk,kk], coef[kk])$test
	    print1 <- rbind(print1, c(NA, NA, NA,
				      temp, x$df[i], 1-pchisq(temp, 1)))
	    }
	else {
	    pname1 <- c(pname1, names(coef)[kk])
	    tempe<- (diag(x$var))[kk]
	    temp <- coef[kk]^2/ tempe
	    print1 <- rbind(print1, cbind(coef[kk], sqrt(tempe),
				      sqrt((diag(x$var2))[kk]), 
				      temp, 1, 1-pchisq(temp, 1)))
	    }
	}

    # Format out the NA's 
    temp <- cbind(format(print1[,1]), format(print1[,2]), 
		       format(print1[,3]),
		       format(round(print1[,4], 2)),
		       format(round(print1[,5], 2)),
		       format(signif(print1[,6], 2)))
    temp <- ifelse(is.na(print1), "", temp)
    dimnames(temp) <- list(substring(pname1,1, maxlabel), 
			     c("coef","se(coef)", "se2", "Chisq","DF","p"))
    prmatrix(temp, quote=FALSE)
	
    #
    # Write out the remaider of the info
    #
    if (nrow(x$var)==length(coef)) 
	    cat("\nScale fixed at",format(x$scale),"\n") 
    else if (length(x$scale)==1) cat ("\nScale=", format(x$scale), "\n")
    else {
	cat("\nScale:\n")
	print(x$scale, ...)
	}

    cat("\nIterations:", x$iter[1], "outer,", x$iter[2], 
	"Newton-Raphson\n")
    if (length(print2)) {
#	cat("Penalized terms:\n")
	for (i in 1:length(print2)) cat("    ", print2[i], "\n")
	}

    logtest <- -2 * (x$loglik[1] - x$loglik[2])
    df <- sum(x$df) - x$idf
#    cat("\n")
    cat("Degrees of freedom for terms=", format(round(x$df,1)), "\n")
#    cat("Loglik (initial,final) = ", format(round(x$loglik,2)),
#	"  Penalty = ", format(x$penalty), "\n")
    cat("Likelihood ratio test=", format(round(logtest, 2)), "  on ",
	round(df,1), " df,", " p=", format(1 - pchisq(logtest, df)),  sep="")

    n <- length(x$linear.predictors)
    omit <- x$na.action
    if (length(omit))
	cat("\n  n=", n, " (", naprint(omit), ")\n", sep="")
    else cat("  n=", n, "\n")
    invisible()
    }
