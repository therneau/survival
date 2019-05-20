# $Id: frailty.brent.S 11166 2008-11-24 22:10:34Z therneau $
#
# Brent's method for finding a maximum
#  If upper and/or lower is given, it transforms x to stay out of trouble
#  during the "bracketing" phase
#
frailty.brent <- function(x, y, lower, upper) {
    n <- length(x)
    if (length(y) != n) stop ("Length mismatch for x and y")
 
    if (n<3) return(mean(x))

    # First, is the solution bracketed?
    # If not, take big steps until it is
    ord <- order(x)
    xx <- x[ord]
    yy <- y[ord]
    best <- (1:n)[yy==max(y)]
    if (length(best) >1) stop("Ties for max(y), I surrender") #fix this later
    if (best==1) {
	new <- xx[1] - 3*(xx[2] - xx[1])
	if (!missing(lower) && !is.null(lower) && new < lower)
		new <- lower + (min(xx[xx>lower])-lower)/10
	return(new)
	}
    if (best==n) {
	new <- xx[n] + 3*(xx[n] - xx[n-1])
	if (!missing(upper) && !is.null(upper) && new > upper)
		new <- upper + (max(xx[xx<upper])-upper)/10
	return(new)
	}

    # Ok, it's bracketed.  Do a quadratic extrapolation
    # Now, these are my best 3 guesses so far
    xx <- xx[(best-1): (best+1)]
    yy <- yy[(best-1): (best+1)]
    temp1 <- (xx[2] -xx[1])^2 *(yy[2]-yy[3]) - (xx[2]-xx[3])^2 * (yy[2]-yy[1])
    temp2 <- (xx[2] -xx[1])   *(yy[2]-yy[3]) - (xx[2]-xx[3])   * (yy[2]-yy[1])
    new <- xx[2] - .5*temp1/temp2

    # if the new guess is outside the bracketing interval, or it is
    #   "bouncing around", then use golden section
    if (new < xx[1] || new > xx[3] ||
	      ( (n>4) && (new-x[n]) > .5*abs(x[n-1]-x[n-2]))) {
        if ((xx[2]-xx[1]) > (xx[3]-xx[2]))  return(xx[2] - .38*(xx[2]-xx[1]))
        else                                return(xx[2] + .32*(xx[3]-xx[2]))
	}
    else return(new)
    }
    
