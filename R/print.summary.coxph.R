print.summary.coxph <-
 function(x, digits = max(getOption('digits')-3, 3),  
             signif.stars = getOption("show.signif.stars"), expand=FALSE, ...) {
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

    if (expand && !is.null(x$cmap)) { # this was a coxphms object
        signif.stars <- FALSE  #work around issue with printCoefmat
        # print it group by group
        tmap <- x$cmap
        cname <- colnames(tmap)
        printed <- rep(FALSE, length(cname))
        for (i in 1:length(cname)) {
            # if multiple colums of tmat are identical, only print that
            #  set of coefficients once
            if (!printed[i]) { # this column hasn't been printed
                j <- apply(tmap, 2, function(x) all(x == tmap[,i])) 
                printed[j] <- TRUE  # mark all that match as 'printed'
               
                tmp2 <- x$coefficients[tmap[,i],, drop=FALSE]
                names(dimnames(tmp2)) <- c(paste(cname[j], collapse=", "), "")
                # restore character row names
                rownames(tmp2) <- rownames(tmap)[tmap[,i]>0]
 
                printCoefmat(tmp2, digits=digits, P.values=TRUE, 
                             has.Pvalue=TRUE, signif.legend=FALSE,
                             signif.stars = signif.stars, ...)

                if (!is.null(x$conf.int)) {
                    tmp2 <- x$conf.int[tmap[,i],, drop=FALSE]
                    rownames(tmp2) <- rownames(tmap)[tmap[,i] >0]
                    names(dimnames(tmp2)) <- c(paste(cname[j], collapse=", "),"")
                    print(tmp2, digits=digits, ...)
                }   
            } 
        }        
        cat("\n States:", paste(paste(seq(along.with =x$states), x$states, sep='= '), 
                               collapse=", "), '\n')
    } else {
        if(!is.null(x$coefficients)) {
            cat("\n")
            printCoefmat(x$coefficients, digits=digits,
                         signif.stars=signif.stars, ...)
        }
        if(!is.null(x$conf.int)) {
            cat("\n")
            print(x$conf.int)
        }
    }       
    cat("\n")

    if (!is.null(x$concordance)) {
        cat("Concordance=", format(round(x$concordance[1],3)),
            " (se =", format(round(x$concordance[2], 3)),")\n")
    }
#    cat("Rsquare=", format(round(x$rsq["rsq"],3)),
#            "  (max possible=", format(round(x$rsq["maxrsq"],3)),
#            ")\n" )

    pdig <- max(1, getOption("digits")-4)  # default it too high IMO
    cat("Likelihood ratio test= ", format(round(x$logtest["test"], 2)), "  on ",
	x$logtest["df"], " df,", "   p=", 
        format.pval(x$logtest["pvalue"], digits=pdig),
        "\n", sep = "")
    cat("Wald test            = ", format(round(x$waldtest["test"], 2)), "  on ",
	x$waldtest["df"], " df,", "   p=", 
        format.pval(x$waldtest["pvalue"], digits=pdig),
	"\n", sep = "")
    cat("Score (logrank) test = ", format(round(x$sctest["test"], 2)), "  on ",
        x$sctest["df"]," df,", "   p=", 
        format.pval(x$sctest["pvalue"], digits=pdig), sep ="")
    if (is.null(x$robscore))
        cat("\n\n")
    else cat(",   Robust = ", format(round(x$robscore["test"], 2)), 
             "  p=", 
             format.pval(x$robscore["pvalue"], digits=pdig), "\n\n", sep="")   

    if (x$used.robust)
	cat("  (Note: the likelihood ratio and score tests",
            "assume independence of\n     observations within a cluster,",
	    "the Wald and robust score tests do not).\n")
    invisible()
}
