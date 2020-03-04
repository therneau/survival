print.coxph <-
 function(x, digits=max(1L, getOption("digits") - 3L), signif.stars=FALSE, ...)
    {
    if (!is.null(cl<- x$call)) {
	cat("Call:\n")
	dput(cl)
	cat("\n")
	}
    if (!is.null(x$fail)) {
	cat(" Coxph failed.", x$fail, "\n")
	return()
	}
    savedig <- options(digits = digits)
    on.exit(options(savedig))

    coef <- x$coefficients
    se <- sqrt(diag(x$var))
    if(is.null(coef) | is.null(se))
        stop("Input is not valid")

    if (is.null(x$naive.var)) {
	tmp <- cbind(coef, exp(coef), se, coef/se,
               pchisq((coef/se)^2, 1, lower.tail=FALSE))
	dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)",
	    "se(coef)", "z", "p"))
	}
    else {
	nse <- sqrt(diag(x$naive.var))
	tmp <- cbind(coef, exp(coef), nse, se, coef/se,
	       pchisq((coef/se)^2, 1, lower.tail=FALSE))
	dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)",
	    "se(coef)", "robust se", "z", "p"))
	}
    
    if (inherits(x, "coxphms")) {
        # print it group by group
        tmap <- x$cmap[-1,,drop=FALSE]  # ignore the intercept (strata)
        cname <- colnames(tmap)
        printed <- rep(FALSE, length(cname))
        for (i in 1:length(cname)) {
            # if multiple colums of tmat are identical, only print that
            #  set of coefficients once
            if (!printed[i]) {
                j <- apply(tmap[-1,, drop=FALSE], 2, 
                           function(x) all(x == tmap[-1,i])) 
                printed[j] <- TRUE

                tmp2 <- tmp[tmap[,i],, drop=FALSE]
                names(dimnames(tmp2)) <- c(paste(cname[j], collapse=", "), "")
                # restore character row names
                rownames(tmp2) <- rownames(tmap)[tmap[,i]>0]
                printCoefmat(tmp2, digits=digits, P.values=TRUE, 
                             has.Pvalue=TRUE,
                             signif.stars = signif.stars, ...)
                cat("\n")
            }       
        }

        cat(" States: ", paste(paste(seq(along=x$states), x$states, sep='= '), 
                               collapse=", "), '\n')
        # cat(" States: ", paste(x$states, collapse=", "), '\n')
        if (FALSE) { # alternate forms, still deciding which I like
            stemp <- x$states
            names(stemp) <- 1:length(stemp)
            print(stemp, quote=FALSE)
        }
    }
    else printCoefmat(tmp, digits=digits, P.values=TRUE, has.Pvalue=TRUE,
                 signif.stars = signif.stars, ...)

    logtest <- -2 * (x$loglik[1] - x$loglik[2])
    if (is.null(x$df)) df <- sum(!is.na(coef))
    else  df <- round(sum(x$df),2)
    cat("\n")
    cat("Likelihood ratio test=", format(round(logtest, 2)), "  on ",
	df, " df,", " p=", 
        format.pval(pchisq(logtest, df, lower.tail=FALSE), digits=digits), 
        "\n",  sep="")
    omit <- x$na.action
    cat("n=", x$n)
    if (!is.null(x$nevent)) cat(", number of events=", x$nevent, "\n")
    else cat("\n")
    if (length(omit))
	cat("\   (", naprint(omit), ")\n", sep="")
    invisible(x)
    }

coef.coxphms <- function(object, type=c("vector", "matrix"), ...) {
    type <- match.arg(type)
    if (type=="matrix") {
        cmap2 <- object$cmap[-1,, drop=FALSE]
        cmat <- 0*cmap2  # all the right names
        cmat[cmap2>0] <- object$coefficient[cmap2]
        attr(cmat, "states") <- object$states
        cmat
    }
    else NextMethod(object, ...)
}

        
