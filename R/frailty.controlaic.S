#  $Id: frailty.controlaic.S 11166 2008-11-24 22:10:34Z therneau $
# Control function to minimize the AIC
#  the optional paramater "caic" chooses corrected aic (default=FALSE)
# n is the "effective" sample size
#

frailty.controlaic <- function(parms, iter, old, n, df, loglik) {
    if (iter==0) {  # initial call
	if (is.null(parms$init)) theta <-0.005
	else theta <- parms$init[1]
	return(list(theta=theta, done=FALSE))
	}
    
    # by default, do the corrected AIC
    if (length(parms$caic)) correct <- parms$caic
    else correct <- FALSE

    if (n < df+2) dfc <- (df -n) + (df+1)*df/2 -1  #avoid pathology
    else          dfc <- -1 + (df+1)/(1- ((df+2)/n))
    if (iter==1) { # Second guess in series
	history <- c(theta=old$theta, loglik=loglik,
		     df=df, aic=loglik-df, aicc=loglik - dfc)
	if (length(parms$init) <2) theta <-1
	else theta <- parms$init[2]
	temp <- list(theta=theta, done=FALSE, history=history)
	return(temp)
	}

    history <- rbind(old$history,c(old$theta, loglik, df, loglik-df, 
				   loglik -dfc))
    if (is.null(parms$trace)) trace <-FALSE
    else trace <- parms$trace
    
    if (iter==2) {  #Third guess
	theta <- mean(history[,1])
	return(list(theta=theta, done=FALSE, history=history))
	}
    #
    # Ok, now we're ready to actually use prior data
    # Now, history has iter rows, each row contains the 
    # value of theta, the Cox PL, the df, aic, and corrected aic
    if (correct) aic <- history[,5]   #use corrected aic for convergence
    else         aic <- history[,4]

    done <- (abs(1- aic[iter]/aic[iter-1]) < parms$eps)
    x <- history[,1]
    
    if (x[iter]== max(aic) && x[iter]==max(x)) 
	    newtheta <- 2* max(x)
    else  newtheta <- frailty.brent(x, aic, lower=parms$lower, 
				    upper=parms$upper)
    
    if (length(parms$trace) && parms$trace) {
	print(history)
	cat("    new theta=", format(newtheta), "\n\n")
	}
    list(theta=newtheta, done=done, history=history)
    }
        
