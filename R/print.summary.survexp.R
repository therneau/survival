print.summary.survexp <- function(x, 
				  digits = max(options()$digits - 4, 3), ...) {
    savedig <- options(digits=digits)
    on.exit(options(savedig))

    if (!is.null(cl<- x$call)) {
	cat("Call: ")
	dput(cl)
	cat("\n")
	}

    omit <- x$na.action
    if (length(omit)) 
	    cat(naprint(omit), "\n")

    mat <- cbind(x$time, x$n.risk,  x$surv)
    cnames <- c("time", "n.risk")

    if (is.matrix(x$surv)) ncurve <- ncol(x$surv)
    else	           ncurve <- 1
    if (ncurve==1) {                 #only 1 curve
	cnames <- c(cnames, "survival")
	if (!is.null(x$std.err)) {
	    if (is.null(x$lower)) {
		mat <- cbind(mat, x$std.err)
		cnames <- c(cnames, "std.err")
	        }
	    else {
		mat <- cbind(mat, x$std.err, x$lower, x$upper)
		cnames <- c(cnames, 'std.err',
			    paste("lower ", x$conf.int*100, "% CI", sep=''),
			    paste("upper ", x$conf.int*100, "% CI", sep=''))
	        }	
	    }
        }
    else cnames <- c(cnames, paste("survival", seq(ncurve), sep=''))

    if (!is.matrix(mat)) mat <- matrix(mat, nrow=1)
    if (!is.null(mat)) {
	dimnames(mat) <- list(NULL, cnames)
	if (is.null(x$strata))
		prmatrix(mat, rowlab=rep("", nrow(mat)))
	else  { #print it out one strata at a time
	    strata <- x$strata
	    for (i in levels(strata)) {
		who <- (strata==i)
		cat("               ", i, "\n")
		if (sum(who) ==1)
			print(mat[who,])
	        else
		    prmatrix(mat[who,], rowlab=rep("", sum(who)))

		cat("\n")
 	        }
	    }
        }
    else 
	stop("There are no observations to print.")
    invisible(x)
    }














