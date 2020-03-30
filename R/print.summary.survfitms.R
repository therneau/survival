print.summary.survfitms <- function(x, 
				  digits = max(options()$digits - 4, 3), ...) {
    savedig <- options(digits=digits)
    on.exit(options(savedig))

    if (!is.null(cl<- x$call)) {
	cat("Call: ")
	dput(cl)
	cat("\n")
	}
    tsum <- function(x) {
        if (is.matrix(x)) rowSums(x) else x
    }
    omit <- x$na.action
    if (length(omit)) 
	    cat(naprint(omit), "\n")

    mat <- cbind(x$time, tsum(x$n.risk), tsum(x$n.event))
    cnames <- c("time", "n.risk", "n.event")

    # If there is only one state, print the estimate, se, and CI at each
    #   time point.  If there are multiples, print just the estimate.
    # If there are multiple strata, or multiple newdata rows (predicted
    #   curves from a multi-state coxph model), then print those sequentially
    #   with a header line between each.
    if (is.null(x$states)) nstate <-1
    else nstate <- length(x$states)

    dd <- dim(x$pstate)
    if (length(dd) ==3 ) {
        if (is.null(x$strata)) group <- rep(paste("data", 1:dd[2]), each=dd[1])
        else group <- c(outer(rep(names(x$strata), x$strata),
                              paste("data", 1:dd[2]), paste, sep=(", ")))
        mat <- mat[rep(1:nrow(mat), dim(x$pstate)[2]), ]
        mat <- cbind(mat, matrix(x$pstate, ncol= dd[3]))
    } else {
        if (is.null(strata)) group <- NULL
        else group <- rep(names(x$strata), x$strata)
        mat <- cbind(mat, x$pstate)
    }

    if (nstate >1) 
        cnames <- c(cnames, paste0("P(", x$states[1:nstate], ")"))
    else {
        cnames <- c(cnames, "P")
	if (!is.null(x$std.err)) {
	    if (is.null(x$lower)) {
		mat <- cbind(mat, as.vector(x$std.err))
		cnames <- c(cnames, "std.err")
            }
	    else {
		mat <- cbind(mat, as.vector(x$std.err), 
                             as.vector(x$lower), as.vector(x$upper))
		cnames <- c(cnames, 'std.err',
			    paste("lower ", x$conf.int*100, "% CI", sep=''),
			    paste("upper ", x$conf.int*100, "% CI", sep=''))
            }	
        }
    }


    if (!is.null(x$start.time)) {
	mat.keep <- mat[,1] >= x$start.time
        if (!any(mat.keep))
            stop(paste("No rows remain using start.time =", x$start.time, "."))
	mat <- mat[mat.keep,,drop=FALSE]
        if (!is.null(group)) group <- group[mat.keep]
    }

    if (nrow(mat) > 0) {
        dimnames(mat) <- list(rep("", nrow(mat)), cnames)
        if (is.null(group)) print(mat)
        else  { #print it out one group at a time
            for (i in unique(group)) {
                who <- (group==i)
                cat("               ", i, "\n")
                print(mat[who,])
                cat("\n")
            }
        }
    } else 
	stop("There are no events to print.  Please use the option ",
	    "censored=TRUE with the summary function to see the censored ",
	    "observations.")
    invisible(x)
    }
