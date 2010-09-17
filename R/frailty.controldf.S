# $Id: frailty.controldf.S 11373 2009-10-28 17:12:59Z therneau $
# A function to calibrate the df
#    very empirical  
# Find the closest 3 points that span the target value
#   We know the function is monotone, so fit the function
#     dy = a * (dx)^p   to the 3 points, where dx and dy are the distance
#     from the leftmost of the three points.
#   This method can fail near a boundary, so use step halving if things don't
#     go well
# On input, parms$df = target degrees of freedom
#           parms$dfs, parms$thetas = known values (usually 0,0)
#           parms$guess = first guess
#
frailty.controldf <- function(parms, iter, old, df) {
    if (iter==0) {  
	theta <- parms$guess
	return(list(theta=theta, done=FALSE, 
		    history=cbind(thetas=parms$thetas, dfs=parms$dfs)))
	}

    eps <- parms$eps
    if (length(eps)==0) eps <- .1

    thetas <- c(old$history[,1], old$theta)
    dfs    <- c(old$history[,2], df)
    nx <- length(thetas)
    if (nx==2) {
	#linear guess based on first two 
	# but try extra hard to bracket the root
	theta <- thetas[1] + (thetas[2]-thetas[1])*(parms$df - dfs[1])/
						    (dfs[2] - dfs[1])
	if (parms$df > df) theta <- theta * 1.5 
	return(list(theta=theta, done=FALSE,
		    history=cbind(thetas=thetas, dfs=dfs), half=0))
	}
    else{
	# Now, thetas= our guesses at theta
	#  dfs = the degrees of freedom for each guess
	done <- (iter>1 &&
		 (abs(dfs[nx]-parms$df) < eps))

	# look for a new minimum
	x <- thetas
	y <- dfs
	target <- parms$df

	# How am I doing
	if ( abs( (y[nx]-target)/(y[nx-1]-target)) > .6) doing.well <- FALSE
	else doing.well <- TRUE
	
	ord <- order(x)
	if ((x[1]-x[2])*(y[1]-y[2]) >0)  y <- y[ord]  #monotone up
	else  { #monotone down
	    y <- -1* y[ord]
	    target <- -target
	    }
	x <- x[ord]

	if (all(y>target)) b1 <- 1     #points 1:3 are the closest then
	else if (all(y<target)) b1 <- nx-2
	else {
	    b1 <- max((1:nx)[y <= target]) #this point below target, next above
	    if (!doing.well && (is.null(old$half) ||  old$half<2)) {
		#try bisection
		if (length(parms$trace) && parms$trace){
		    print(cbind(thetas=thetas, dfs=dfs))
		    cat("  bisect:new theta=" , format( mean(x[b1+0:1])), 
			"\n\n")
		    }
		return(list(theta= mean(x[b1+0:1]),done=done, 
			      history=cbind(thetas=thetas, dfs=dfs), 
				            half=max(old$half, 0) +1))
		}
	    # use either b1,b1+1,b1+2 or  b1-1, b1, b1+1, whichever is better
	    #  better = midpoint of interval close to the target

	    if ((b1+1)==nx ||
		(b1>1 &&  ((target -y[b1]) < (y[b1+1] -target))))
		    b1 <- b1-1
	    }

	#now have the best 3 points
	# fit them with a power curve anchored at the leftmost one
	b2 <- b1 + 1:2	
	xx <- log(x[b2] - x[b1])
	yy <- log(y[b2] - y[b1])
	power <- diff(yy)/diff(xx)
	a <- yy[1] - power*xx[1]
	newx <- (log(target -y[b1]) - a)/power
	if (length(parms$trace) && parms$trace){
	    print(cbind(thetas=thetas, dfs=dfs))
	    cat("  new theta=" , format(x[b1] + exp(newx)), "\n\n")
	    }
	list(theta=x[b1] + exp(newx), done=done, 
	     history=cbind(thetas=thetas, dfs=dfs), half=0)
	}
    }

