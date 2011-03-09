# $Id $
print.summary.coxph <-
 function(x, digits = max(getOption('digits')-3, 3),  
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
    if (!is.null(x$nevent)) cat(", number of events=", x$nevent, "\n")
    else cat("\n")
    if (length(omit))
	cat("   (", naprint(omit), ")\n", sep="")

    if (nrow(x$coef)==0) {   # Null model
	cat ("   Null model\n")
	return()
        }


    if(!is.null(x$coefficients)) {
        cat("\n")
        if (is.R()) printCoefmat(x$coefficients, digits=digits,
                                 signif.stars=signif.stars, ...)
        else prmatrix(x$coefficients)
        }
    if(!is.null(x$conf.int)) {
        cat("\n")
        prmatrix(x$conf.int)
        }
    cat("\n")

    if (!is.null(x$concordance)) {
        cat("Concordance=", format(round(x$concordance[1],3)),
            " (se =", format(round(x$concordance[2], 3)),")\n")
    }
    cat("Rsquare=", format(round(x$rsq["rsq"],3)),
            "  (max possible=", format(round(x$rsq["maxrsq"],3)),
            ")\n" )

    cat("Likelihood ratio test= ", format(round(x$logtest["test"], 2)), "  on ",
	x$logtest["df"], " df,", "   p=", format(x$logtest["pvalue"]),
        "\n", sep = "")
    cat("Wald test            = ", format(round(x$waldtest["test"], 2)), "  on ",
	x$waldtest["df"], " df,", "   p=", format(x$waldtest["pvalue"]),
	"\n", sep = "")
    cat("Score (logrank) test = ", format(round(x$sctest["test"], 2)), "  on ",
        x$sctest["df"]," df,", "   p=", format(x$sctest["pvalue"]), sep ="")
    if (is.null(x$robscore))
        cat("\n\n")
    else cat(",   Robust = ", format(round(x$robscore["test"], 2)), 
             "  p=", format(x$robscore["pvalue"]), "\n\n", sep="")   

    if (x$used.robust)
	cat("  (Note: the likelihood ratio and score tests",
            "assume independence of\n     observations within a cluster,",
	    "the Wald and robust score tests do not).\n")
    invisible()
}
