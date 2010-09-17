# $Id: frailty.controlgauss.S 11166 2008-11-24 22:10:34Z therneau $
#
# The control function for REML on a gaussian
#
frailty.controlgauss <- function(opt, iter, old, fcoef, trH, loglik){
    if (iter==0) {
	# initial call
	#  Because of how the iteration works, 0 is not a useful trial value
	if (!is.null(opt$theta)) theta <- opt$theta  #fixed theta case
	else {
	    if (is.null(opt$init)) theta <- 1
	    else theta <- opt$init[1]
	    }
	list(theta=theta)
        }
    
    else {
	if (is.null(opt$trace)) trace <-FALSE
	else trace <- opt$trace

	nfrail <- length(fcoef)
	fsum   <- sum(fcoef^2)
	theta <- old$theta
	resid <- fsum/(nfrail - trH/theta) - theta

	# save history of the iteration, and get the next theta
	if (iter==1) {
	    history <- c(theta=theta, resid=resid, fsum=fsum, trace=trH)
	    if (is.null(opt$init )) {
		if (resid>0)  theta <- theta*3
		else          theta <- theta/3  	        }
	    else theta <- opt$init[2]
	    list(theta=theta, done=FALSE, history=history)
	    }
	else {
            history <- rbind(old$history,
			     as.vector(c(theta, resid, fsum, trH)))
	    if (iter ==2) {
		if (all(history[,2] > 0))     theta <- history[2,1]*2
		else if (all(history[,2] <0)) theta <- history[2,1]/2
		else    		theta <- mean(history[1:2,1])
                if (trace) {
		    print(history)
		    cat("    new theta=", theta, "\n\n")
		    }
		list(theta=theta, done=FALSE, history=history)
		}
	    else {
		done <- (abs(history[iter,2]) < opt$eps)
		ord <- order(history[,1])
		tempy <- history[ord,2]  #x & y from left to right
		tempx <- history[ord,1]
		# make sure we have one positve and one negative y value
		#  y must be positive near 0, and negative for large x
		if (all(tempy>0))  newtheta <- 2*max(tempx)
		else if (all(tempy<0)) newtheta <- .5 * min(tempx)
		else{
		    #find the latest point, and one on each side of 0
		    b1 <- (1:iter)[ord==iter]
		    if (b1==1) b1 <-2
		    else if (b1==iter) b1 <- iter-1

		    # Brent's formula, straight from Numerical Recipies
		    # Why, you may ask, don't we use the uniroot() function
		    #  which is built into S, and implements Brent's method?
		    # Because all we want is the next guess for x.  The interal
		    #  loop of coxph is calling us, not the other way around.
		    guess <- history[iter- (2:0),1]
		    R <- tempy[b1]/ tempy[b1+1]
		    S <- tempy[b1]/ tempy[b1-1]
		    U <- R/S
		    P <- S* (U*(R-U)*(tempx[b1+1]-tempx[b1]) - 
			    (1-R)*(tempx[b1]-tempx[b1-1]))
		    Q <- (U-1)*(R-1)*(S-1)
		    newtheta <- tempx[b1] + P/Q
		    # if the new guess is outside the brackets, do a binomial
		    #   search step
		    if (newtheta > tempx[b1+1]) newtheta <- mean(tempx[b1+0:1])
		    if( newtheta < tempx[b1-1]) newtheta <- mean(tempx[b1-0:1])
		    }
                if (trace) {
		    print(history)
		    cat("    new theta=", format(newtheta), "\n\n")
		    }
		list(theta=newtheta, done=done, history=history)
		}
	    }
        }
    }
	
