# $Id: lines.aareg.S 11166 2008-11-24 22:10:34Z therneau $
lines.aareg <- function(x, se=FALSE, maxtime, type='s', ...) {
    if (!inherits(x, 'aareg')) stop ("Must be an aareg object")

    if (missing(maxtime)) keep <- 1:length(x$time)
    else		  keep <- 1:sum(x$time <= maxtime)

    if (is.matrix(x$coefficient) && ncol(x$coefficient)>1) {
	yy <- apply(x$coefficient[keep,], 2,cumsum)
	yy <- rbind(0,yy)  # make the plot start at 0,0
	if (se) {
	    if (!is.null(x$dfbeta)) {
		# There was a cluster term, so use the robust variance
		#  dfbeta will be of dimension (n, nvar, n-unique-times)
		# The first variance increment is apply(dfbeta[,,1]^2,2,sum)
		#              second is          apply(dfbeta[,,2]^2,2,sum)
                #              ... , apply(dfbeta[,,ndeath].....
		# By being sneaky, it can be done quickly
		dd <- dim(x$dfbeta)
		keep2 <- 1:length(unique(x$time[keep]))
		temp <- matrix(x$dfbeta[,,keep2], nrow=dd[1])
		se.increment <- matrix(apply(temp^2, 2, sum), nrow=dd[2])
		se.yy <- sqrt(apply(t(se.increment), 2, cumsum))
		}
	    else se.yy <- sqrt(apply(x$coefficient[keep,]^2, 2,cumsum)) 
	    se.yy <- rbind(0, se.yy)
	    }
	ncurve <- ncol(yy)
	}
    
    else {
	# this is the branch most often called, when someone has done
	#   plot(fit[3]), so that only 1 coefficient remains
	yy <- cumsum(c(0, x$coefficient[keep]))
	if (se) {
	    if (!is.null(x$dfbeta)) {
		dd <- dim(x$dfbeta)
		keep2 <- 1:length(unique(x$time[keep]))
		temp <- matrix(x$dfbeta[,,keep2], nrow=dd[1])
		se.yy <- sqrt(cumsum(c(0, apply(temp^2, 2, sum))))
		}
	    else se.yy <- sqrt(cumsum(c(0, x$coefficient[keep]^2)))
	    }
	ncurve <- 1
	}
   
    xx <- c(0, x$time[keep])
	
    # There may be multiplicities in x$times.  Only plot the last of
    #  each of them
    indx <- 1 + length(xx) - rev(match(unique(rev(xx)), rev(xx)))
    xx <- xx[indx]
    yy <- as.matrix(yy)[indx,]

    if (se) { 
	if (is.null(x$dfbeta)) se.yy<- as.matrix(se.yy)[indx,]
	yy <- cbind(yy, yy + 1.96*se.yy,
			 yy - 1.96*se.yy)
	if (ncurve >1) {
	    for (i in 1:ncurve) {
		j <- c(i, i+ncurve, i+2*ncurve)
		matlines(xx, yy[,j], type=type, ..., col=1, lty=c(1,2,2))
		}
	    }
	else matlines(xx, yy, type=type, ...,  col=1, lty=c(1,2,2),)
	}
    else {
	matlines(xx, yy, type=type, ...,  xlab='Time')
	}
    }
