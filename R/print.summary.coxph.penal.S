print.summary.coxph.penal <-
 function(x,  digits = max(options()$digits - 3, 3),
           signif.stars = getOption("show.signif.stars"), ...) {
    if (!is.null(x$call)) {
	cat("Call:\n")
	dput(x$call)
	cat("\n")
	}
    if (!is.null(x$fail)) {
	cat(" Coxreg failed.", x$fail, "\n")
	return()
	}
    savedig <- options(digits = digits)
    on.exit(options(savedig))

    omit <- x$na.action
    cat("  n=", x$n)
    if (!is.null(x$nevent)) 
        cat(", number of events=", x$nevent, "\n")
    else cat("\n")
    if (length(omit))
	cat("   (", naprint(omit), ")\n\n", sep="")
    else cat("\n")

    # Format out the NA in the coef matrix
    print1 <- x$coefficients
    temp <- cbind(format(print1[,1]), format(print1[,2]), 
		       format(print1[,3]),
		       format(round(print1[,4], 2)),
		       format(round(print1[,5], 2)),
		       format(signif(print1[,6], 2)))
    temp <- ifelse(is.na(print1), "", temp)
    dimnames(temp) <- dimnames(print1)
    print(temp, quote=FALSE)

    if(length(x$conf.int) >0 ) {
        cat("\n")
        print(x$conf.int)
        }
    logtest <- -2 * (x$loglik[1] - x$loglik[2])
    sctest <- x$score

    cat("\nIterations:", x$iter[1], "outer,", x$iter[2], 
        "Newton-Raphson\n")
    if (length(x$print2)) {
        for (i in 1:length(x$print2)) cat("    ", x$print2[i], "\n")
        }
    if (is.null(x$df)) df <- sum(!is.na(coef))
    else  df <- round(sum(x$df),2)
    cat("Degrees of freedom for terms=", format(round(x$df,1)), "\n")
    if (!is.null(x$concordance)) {
        cat("Concordance=", format(round(x$concordance[1],3)),
            " (se =", format(round(x$concordance[2], 3)),")\n")
    }
    cat("Likelihood ratio test= ", format(round(logtest, 2)), "  on ",
	df, " df,", "   p=", format(1 - pchisq(logtest, df)),
	"\n", sep = "")
    if (!is.null(x$wald.test))
        cat("Wald test            = ", format(round(x$wald.test, 2)), 
	    "  on ", df, " df,   p=",
	    format(1 - pchisq(x$wald.test, df)), sep = "")
    if (!is.null(x$score))
        cat("\nScore (logrank) test = ", format(round(sctest, 2)), "  on ", df,
            " df,", "   p=", format(1 - pchisq(sctest, df)), sep ="") 
    if (is.null(x$rscore)) cat("\n")
    else cat(",   Robust = ", format(round(x$rscore, 2)), 
	   "  p=", format(1 - pchisq(x$rscore, df)), "\n", sep="")   

    invisible()
    }
