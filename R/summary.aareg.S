# The summary routine for aareg models.
# A lot of the work below relates to one particular issue: the coeffients
#  of an aareg model often get "wild" near the end (at the largest times).
#  So, a common case is to 
#       fit the model (very slow)
#       look at the printout -- Hmmm x1 is significant, x2 not, ...., why?
#       look at plot(fit) and
#           oh my gosh, I should have cut the time scale off at 520 days
#
# This routine allows one to do that.  If maxtime is given, the overall
#    test statistic is re-computed.  One consequence is that lots of the
#    intermediate material from the fit had to be included in the aareg
#    object.
# The "variance" based weighting for a test is not allowed, because it would
#    have meant an awful lot more stuff to pass, lots more work, for a test
#    that is rarely used.
#
summary.aareg <- function(object, maxtime, test=c('aalen', 'nrisk'), scale=1,...) {
    if (!inherits(object, 'aareg')) stop ("Must be an aareg object")

    if (missing(test)) test <- object$test
    test <- match.arg(test)

    if (!missing(maxtime)) ntime <- sum(object$time <= maxtime)
    else 		   ntime <- nrow(object$coefficient)
    times <- object$time[1:ntime]

    if (test=='aalen') {
        twt <- (as.matrix(object$tweight))[1:ntime,]
        scale <- apply(twt, 2, sum)/scale
        }
    else {
        twt <-  object$nrisk[1:ntime]
        scale <- ntime/scale
        }

    # Compute a "slope" for each line, using appropriate weighting
    #  Since this is a single variable model, no intercept, I
    #  don't need to call lm.wfit!
    tx <- as.matrix(twt * object$coefficient[1:ntime,])
    ctx  <- apply(tx, 2, cumsum)
    if (is.matrix(twt) && ncol(twt) >1) tempwt <- apply(twt*times^2, 2, sum)
    else                                tempwt <- sum(twt*times^2)
    if (ncol(ctx) >1)
            slope<- apply(ctx* times, 2, sum)/ tempwt
    else    slope <- sum(ctx*times) / tempwt

    if (!missing(maxtime) || object$test != test) {
	# Compute the test statistic
	test.stat <- apply(tx, 2, sum)  #sum of tested coefficients
	test.var  <- t(tx) %*% tx	#std Poisson, ind coefficients variance
	
	if (!is.null(object$dfbeta)) {
	    dd <- dim(object$dfbeta)
	    indx <- match(unique(times), object$times)
	    influence <- matrix(0, dd[1], dd[2])
	    for (i in 1:length(indx)) {
		if (test=='aalen') 
		     influence <- influence + object$dfbeta[,,i] %*% 
			                            diag(twt[indx[i],])
		else influence <- influence + object$dfbeta[,,i]* object$nrisk[indx[i]]
		}
	    if (!is.null(object$cluster)) influence <- rowsum(influence, cluster)
	    test.var2 <- t(influence) %*% influence
	    }
        else test.var2 <- NULL
	}
    else {  #use the value that was passed in
	test.stat <- object$test.statistic
	test.var  <- object$test.var
	test.var2 <- object$test.var2   #NULL if dfbeta option was false
	}

    # create the matrix for printing out
    #  The chisquare test does not include the intercept
    se1 <- sqrt(diag(test.var))
    if (is.null(test.var2)) {
	mat <- cbind(slope, test.stat/scale, se1/scale, test.stat/se1, 
		     2*pnorm(-abs(test.stat/se1)))
	dimnames(mat) <- list((dimnames(object$coefficient)[[2]]),
			      c("slope", "coef", "se(coef)", "z", "p"))
	chi <- test.stat[-1] %*% solve(test.var[-1,-1],test.stat[-1]) 
	}
    else {
	se2 <- sqrt(diag(test.var2))
	mat <- cbind(slope, test.stat/scale, se1/scale, se2/scale, 
                     test.stat/se2, 2*pnorm(-abs(test.stat/se2)))
	dimnames(mat) <- list((dimnames(object$coefficient)[[2]]),
			 c("slope", "coef", "se(coef)", "robust se", "z", "p"))
	chi <- test.stat[-1] %*% solve(test.var2[-1,-1], test.stat[-1]) 
	}

    temp <- list(table=mat, test=test, test.statistic=test.stat, 
		 test.var=test.var, test.var2=test.var2, chisq=chi,
		 n = c(object$n[1], length(unique(times)), object$n[3]))

    
    if (is.R()) class(temp) <- 'summary.aareg'
    else      oldClass(temp) <- 'summary.aareg'
    temp
    }

print.summary.aareg <- function(x, ...) {
    print(signif(x$table,3))
    chi <- x$chisq
    df  <- length(x$test.statistic) -1
    cat("\nChisq=", format(round(chi,2)), " on ", df, " df, p=",
	            signif(1- pchisq(chi, df),2), 
	            "; test weights=", x$test, "\n", sep='')
    invisible(x$table)
    }
