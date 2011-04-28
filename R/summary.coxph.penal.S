# $Id $
summary.coxph.penal <-
 function(object, conf.int = 0.95, scale = 1, terms=FALSE, maxlabel=25,
			digits = max(options()$digits - 4, 3),...) {
    if (!is.null(object$call)) {
	cat("Call:\n")
	dput(object$call)
	cat("\n")
	}
    if (!is.null(object$fail)) {
	cat(" Coxreg failed.", object$fail, "\n")
	return()
	}
    savedig <- options(digits = digits)
    on.exit(options(savedig))

    omit <- object$na.action
    omit <- object$na.action
    cat("  n=", object$n)
    if (!is.null(object$nevent)) 
        cat(", number of events=", object$nevent, "\n")
    else cat("\n")
    if (length(omit))
	cat("   (", naprint(omit), ")\n\n", sep="")
    else cat("\n")

    coef <- object$coefficients
    if (length(coef)==0 && length(object$frail)==0)
            stop("Penalized summary function can't be used for a null model")

    if (length(coef) > 0) {
	nacoef <- !(is.na(coef))          #non-missing coefs
	coef2 <- coef[nacoef]
	if(is.null(coef) | is.null(object$var))
		stop("Input is not valid")
	se <- sqrt(diag(object$var))
	}
    #
    # Map terms to special print functions, and the list of iteration histories
    #
    pterms <- object$pterms
    nterms <- length(pterms)
    npenal <- sum(pterms>0)
    print.map <- rep(0,nterms)
    if (!is.null(object$printfun)) {
	temp <- unlist(lapply(object$printfun, is.null))  #which ones are missing
	print.map[pterms>0] <- (1:npenal) * (!temp)
	}

    # Tedious, but build up the coef matrix a term at a time
    print1 <- NULL
    pname1 <- NULL
    if (is.null(object$assign2)) alist <- object$assign[-1]
    else alist <- object$assign2

    print2 <- NULL
    for (i in 1:nterms) {
	kk <- alist[[i]]
	if (print.map[i] >0) {
	    j <- print.map[i]	
	    if (pterms[i]==2) 
		 temp <- (object$printfun[[j]])(object$frail, object$fvar, ,
				  object$df[i], object$history[[j]])
	    else temp <- (object$printfun[[j]])(coef[kk], object$var[kk,kk], 
					   object$var2[kk,kk], 
					   object$df[i], object$history[[j]])
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
	    temp <- coxph.wtest(object$var[kk,kk], coef[kk])$test
	    print1 <- rbind(print1, c(NA, NA, NA,
				      temp, object$df[i], 1-pchisq(temp, 1)))
	    }
	else {
	    pname1 <- c(pname1, names(coef)[kk])
	    tempe<- (diag(object$var))[kk]
	    temp <- coef[kk]^2/ tempe
	    print1 <- rbind(print1, cbind(coef[kk], sqrt(tempe),
				      sqrt((diag(object$var2))[kk]), 
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

    if(conf.int & length(coef) >0 ) {
        z <- qnorm((1 + conf.int)/2, 0, 1)
        coef <- coef * scale
        se <- se * scale
        tmp <- cbind(exp(coef), exp(-coef), exp(coef - z * se),
            exp(coef + z * se))
        dimnames(tmp) <- list(substring(names(coef),1, maxlabel), 
			      c("exp(coef)", "exp(-coef)",
            paste("lower .", round(100 * conf.int, 2), sep = ""),
            paste("upper .", round(100 * conf.int, 2), sep = "")))
        cat("\n")
        prmatrix(tmp)
        }
    logtest <- -2 * (object$loglik[1] - object$loglik[2])
    sctest <- object$score

    cat("\nIterations:", object$iter[1], "outer,", object$iter[2], 
        "Newton-Raphson\n")
    if (length(print2)) {
        for (i in 1:length(print2)) cat("    ", print2[i], "\n")
        }
    if (is.null(object$df)) df <- sum(!is.na(coef))
    else  df <- round(sum(object$df),2)
    cat("Degrees of freedom for terms=", format(round(object$df,1)), "\n")
    if (!is.null(object$concordance)) {
        if (is.matrix(object$concordance)) temp <- colSums(object$concordance)
        else temp <- object$concordance
        temp2 <- c("concordance"= (temp[1] + temp[3]/2)/
                              sum(temp[1:3]), "se"= temp[5]/(2*sum(temp[1:3])))
        cat("Concordance=", format(round(temp2[1],3)),
            " (se =", format(round(temp2[2], 3)),")\n")
    }
    cat("Rsquare=", format(round(1-exp(-logtest/object$n),3)),
	"  (max possible=", format(round(1-exp(2*object$loglik[1]/object$n),3)),
	")\n" )
    cat("Likelihood ratio test= ", format(round(logtest, 2)), "  on ",
	df, " df,", "   p=", format(1 - pchisq(logtest, df)),
	"\n", sep = "")
    if (!is.null(object$wald.test))
        cat("Wald test            = ", format(round(object$wald.test, 2)), 
	    "  on ", df, " df,   p=",
	    format(1 - pchisq(object$wald.test, df)), sep = "")
    if (!is.null(object$score))
        cat("\nScore (logrank) test = ", format(round(sctest, 2)), "  on ", df,
            " df,", "   p=", format(1 - pchisq(sctest, df)), sep ="") 
    if (is.null(object$rscore)) cat("\n")
    else cat(",   Robust = ", format(round(object$rscore, 2)), 
	   "  p=", format(1 - pchisq(object$rscore, df)), "\n", sep="")   

    invisible()
    }
