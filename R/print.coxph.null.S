# $Id: print.coxph.null.S 11166 2008-11-24 22:10:34Z therneau $
print.coxph.null <-
 function(x, digits=max(options()$digits - 4, 3), ...)
    {
    if (!is.null(cl<- x$call)) {
	cat("Call:  ")
	dput(cl)
	cat("\n")
	}

    cat("Null model\n  log likelihood=", format(x$loglik), "\n")
    omit <- x$na.action
    if (length(omit))
	cat("  n=", x$n, " (", naprint(omit), ")\n",
				sep="")
    else cat("  n=", x$n, "\n")
    }
