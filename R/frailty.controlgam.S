# $Id: frailty.controlgam.S 11166 2008-11-24 22:10:34Z therneau $
#
# The control function for a single Gamma frailty term.
#
frailty.controlgam <- function(opt, iter, old, group, status, loglik){
    
    if (iter==0) {
	# initial call
	if (!is.null(opt$theta)) theta <- opt$theta  #fixed theta case
	else {
	    if (is.null(opt$init)) theta <- 0   #no initial value -- use 0
	    else theta <- opt$init[1]
	    }
	list(theta=theta)
        }
    
    else {
	if (is.null(opt$trace)) trace <-FALSE
	else trace <- opt$trace

	theta <- old$theta
	
	#compute correction to the loglik
	if (theta==0) correct <- 0
	else {
	    if (is.matrix(group)) group <- c(group %*% 1:ncol(group))
	    d <- tapply(status,group,sum)
	    correct <- frailty.gammacon(d, 1/theta)
	    }

	if (!is.null(opt$theta)) # fixed theta case
	    list(theta=theta, done=TRUE, c.loglik=loglik + correct)
	else {
	    # save history of the iteration, and get the next theta
	    if (iter==1) history <- c(theta=theta, loglik=loglik, 
		  c.loglik=loglik + correct)
	    else history <- rbind(old$history, 
			 as.vector(c(theta, loglik, 
				     loglik + correct)))
	
	    if (iter==1) {
		if (is.null(opt$init )) theta <-1
		else                    theta <- opt$init[2]
		list(theta=theta, done=FALSE, history=history,
		     c.loglik= loglik+correct)
	        }
	    else if (iter ==2) {
		if (history[2,3] < (history[1,3] +1)) 
			theta <- mean(history[1:2,1])
		else    theta <- 2*history[2,1]

                if (trace) {
		    print(history)
		    cat("    new theta=", theta, "\n\n")
		    }
		list(theta=theta, done=FALSE, history=history,
		     c.loglik= loglik+correct)
		}
	    else {
		#Now, history has iter rows, each row contains the value
	        # of theta, the Cox PL, and the full LL
	        done <- (abs(1- history[iter,3]/history[iter-1,3]) < opt$eps)
		x <- history[,1]
		y <- history[,3]

		if (y[iter]== max(y) && x[iter]==max(x)) newtheta <- 2* max(x)
		else  newtheta <- frailty.brent(sqrt(x), y, lower=0)^2
                if (trace) {
		    print(history)
		    cat("    new theta=", format(newtheta), "\n\n")
		    }
		list(theta=newtheta, done=done, history=history, 
		     c.loglik = loglik + correct)
		}
	    }
        }
    }
	
