# $Id: plot.survfit.S 11377 2009-12-14 22:59:56Z therneau $
plot.survfit<- function(x, conf.int,  mark.time=TRUE,
			mark=3,col=1,lty=1, lwd=1, cex=1, log=FALSE,
			xscale=1, yscale=1, 
			firstx=0, firsty=1,
			xmax, ymin=0,
			fun, xlab="", ylab="", xaxs='S', ...) {

    dotnames <- names(list(...))
    if (any(dotnames=='type'))
        stop("The graphical argument 'type' is not allowed")
    
    if (is.logical(log)) {
	logy <- log
	logx <- FALSE
	if (logy) logax <- 'y'
	else      logax <- ""
        }
    else {
	logy <- (log=='y' || log=='xy')
	logx <- (log=='x' || log=='xy')
	logax  <- log
        }

    if (missing(firstx)) {
	if (!is.null(x$start.time)) 
	     firstx <- x$start.time
	else {
            if (logx || (!missing(fun) && is.character(fun) && fun=='cloglog'))
                firstx <- min(x$time[x$time>0])
            else      firstx <- min(0, x$time)
            }
	}
    firstx <- firstx/xscale

    # The special x axis style only applies when firstx is not given
    if (missing(xaxs) && firstx!=0) xaxs<- par('xaxs')  # use the default

    if (!inherits(x, 'survfit'))
	    stop("First arg must be the result of survfit")

    if (missing(conf.int)) {
	if (is.null(x$strata) && !is.matrix(x$surv)) conf.int <-TRUE
	else conf.int <- FALSE
        }

    if (is.null(x$strata)) {
	nstrat <- 1
	stemp <- rep(1, length(x$time))
        }
    else {
	nstrat <- length(x$strata)
	stemp <- rep(1:nstrat, x$strata)
        }

    ssurv <- x$surv
    stime <- x$time
    supper <- x$upper
    slower <- x$lower
    if (!missing(xmax) && any(x$time>xmax)) {
	# prune back the survival curves
	# I need to replace x's over the limit with xmax, and y's over the
	#  limit with either the prior y value or firsty
	keepx <- keepy <- NULL  # lines to keep
	yzero <- NULL           # if all points on a curve are < xmax
	tempn <- table(stemp)
	offset <- cumsum(c(0, tempn))
	for (i in 1:nstrat) {
	    ttime <-stime[stemp==i]
	    if (all(ttime <= xmax)) {
		keepx <- c(keepx, 1:tempn[i] + offset[i])
		keepy <- c(keepy, 1:tempn[i] + offset[i])
	        }
	    else {
		bad <- min((1:tempn[i])[ttime>xmax])
		if (bad==1)  {
		    keepy <- c(keepy, 1+offset[i])
		    yzero <- c(yzero, 1+offset[i])
		    } 
		else  keepy<- c(keepy, c(1:(bad-1), bad-1) + offset[i])
		keepx <- c(keepx, (1:bad)+offset[i])
		stime[bad+offset[i]] <- xmax
		x$n.event[bad+offset[i]] <- 1   #don't plot a tick mark
	        }
	    }	

	# ok, now actually prune it
	stime <- stime[keepx]
	stemp <- stemp[keepx]
	x$n.event <- x$n.event[keepx]
	if (is.matrix(ssurv)) {
	    if (length(yzero))
		    ssurv[yzero,] <- firsty
	    ssurv <- ssurv[keepy,,drop=FALSE]
	    if (!is.null(supper)) {
		if (length(yzero)) supper[yzero,] <- slower[yzero,] <- firsty
		supper <- supper[keepy,,drop=FALSE]
		slower <- slower[keepy,,drop=FALSE]
	        }
	    }
	else {
	    if (length(yzero)) ssurv[yzero] <- firsty
	    ssurv <- ssurv[keepy]
	    if (!is.null(supper)) {
		if (length(yzero)) supper[yzero] <- slower[yzero] <- firsty
		supper <- supper[keepy]
		slower <- slower[keepy]
  	        }
	    }
        }
    stime <- stime/xscale
    	
    if (!missing(fun)) {
	if (is.character(fun)) {
	    tfun <- switch(fun,
			   'log' = function(x) x,
			   'event'=function(x) 1-x,
			   'cumhaz'=function(x) -log(x),
			   'cloglog'=function(x) log(-log(x)),
			   'pct' = function(x) x*100,
			   'logpct'= function(x) 100*x,
			   stop("Unrecognized function argument")
			   )
	    if (fun=='log'|| fun=='logpct') logy <- TRUE

	    if (fun=='cloglog') {
		logx <- TRUE
		if (logy) logax <- 'xy'
		else logax <- 'x'
	        }
	    }
	else if (is.function(fun)) tfun <- fun
        else stop("Invalid 'fun' argument")
	
	ssurv <- tfun(ssurv )
	if (!is.null(supper)) {
	    supper <- tfun(supper)
	    slower <- tfun(slower)
	    }
	firsty <- tfun(firsty)
	ymin <- tfun(ymin)
        }

    if (is.null(x$n.event)) mark.time <- FALSE   #expected survival curve

    # set default values for missing parameters
    if (is.matrix(ssurv)) ncurve <- nstrat * ncol(ssurv)
    else 		  ncurve <- nstrat

    mark <- rep(mark, length.out=ncurve)
    col  <- rep(col, length.out=ncurve)
    lty  <- rep(lty, length.out=ncurve)
    lwd  <- rep(lwd, length.out=ncurve)

    if (is.numeric(mark.time)) mark.time <- sort(mark.time)

    # Do axis range computations
    if (xaxs=='S') {
	#special x- axis style for survival curves
	xaxs <- 'i'  #what S thinks
	tempx <- max(stime) * 1.04
        }
    else tempx <- max(stime)
    tempx <- c(firstx, tempx, firstx)

    if (logy) {
	tempy <-  range(ssurv[is.finite(ssurv)& ssurv>0])
	if (tempy[2]==1) tempy[2] <- .99
	if (any(ssurv==0)) {
	    tempy[1] <- tempy[1]*.8
	    ssurv[ssurv==0] <- tempy[1]
	    if (!is.null(supper)) {
		supper[supper==0] <- tempy[1]
		slower[slower==0] <- tempy[1]
	        }
	    }
	tempy <- c(tempy, firsty)
        }
    else tempy <- c(range(ssurv[is.finite(ssurv)] ), firsty)
    
    if (missing(fun)) {
	tempx <- c(tempx, firstx)
	tempy <- c(tempy, ymin)
        }
    #
    # Draw the basic box
    #
    plot(tempx, tempy*yscale, type='n', log=logax,
	                  xlab=xlab, ylab=ylab, xaxs=xaxs,...)

    if(yscale != 1) {
	if (logy) par(usr =par("usr") -c(0, 0, log10(yscale), log10(yscale))) 
	else par(usr =par("usr")/c(1, 1, yscale, yscale))   
        }

    #
    # put up the curves one by one
    #   survfit has already put them into the "right" order
    dostep <- function(x,y) {
	if (is.na(x[1] + y[1])) {
	    x <- x[-1]
	    y <- y[-1]
	    }
	n <- length(x)
	if (n > 2) {
	    # replace verbose horizonal sequences like
	    # (1, .2), (1.4, .2), (1.8, .2), (2.3, .2), (2.9, .2), (3, .1)
            # with (1, .2), (3, .1).  They are slow, and can smear the looks
	    # of the line type.
	    dupy <- c(!duplicated(y)[-n], TRUE)
	    n2 <- sum(dupy)

	    #create a step function
	    xrep <- rep(x[dupy], c(1, rep(2, n2-1)))
	    yrep <- rep(y[dupy], c(rep(2, n2-1), 1))
	    list(x=xrep, y=yrep)
            }
	else if (n==1) list(x=x, y=y)
	else           list(x=x[c(1,2,2)], y=y[c(1,1,2)])
        }

    i <- 0
    xend <- NULL
    yend <- NULL

    for (j in unique(stemp)) {
	who <- (stemp==j)
	xx <- c(firstx, stime[who])
	nn <- length(xx)
	if (x$type == 'counting') {
	    deaths <- c(-1, x$n.censor[who])
	    zero.one <- 1
	    }
	else if (x$type == 'right') {
	    deaths <- c(-1, x$n.event[who])
	    zero.one <- 0
	    }
	if (is.matrix(ssurv)) {
	    for (k in 1:ncol(ssurv)) {
		i <- i+1
		yy <- c(firsty, ssurv[who,k])
		lines(dostep(xx, yy), lty=lty[i], col=col[i], lwd=lwd[i]) 

		if (is.numeric(mark.time)) {
		    indx <- mark.time
		    for (k in seq(along.with=mark.time))
			indx[k] <- sum(mark.time[k] > xx)
		    points(mark.time[indx<nn], yy[indx[indx<nn]],
			   pch=mark[i],col=col[i],cex=cex)
		    }
		else if (mark.time && any(deaths==zero.one)) {
		    points(xx[deaths==zero.one], 
			   yy[deaths==zero.one],
			   pch=mark[i],col=col[i],cex=cex)
		    }

		xend <- c(xend,max(xx))
		yend <- c(yend,min(yy))

		if (conf.int && !is.null(supper)) {
		    if (ncurve==1) lty[i] <- lty[i] +1
		    yy <- c(firsty, supper[who,k])
		    lines(dostep(xx,yy), lty=lty[i], col=col[i], lwd=lwd[i])
		    yy <- c(firsty, slower[who,k])
		    lines(dostep(xx,yy), lty=lty[i], col=col[i], lwd=lwd[i])
		    }
	        }
	    }
	else {
	    i <- i+1
	    yy <- c(firsty, ssurv[who])
	    lines(dostep(xx, yy), lty=lty[i], col=col[i], lwd=lwd[i])

	    if (is.numeric(mark.time)) {
		indx <- mark.time
		for (k in seq(along=mark.time))
		    indx[k] <- sum(mark.time[k] > xx)
		points(mark.time[indx<nn], yy[indx[indx<nn]],
		       pch=mark[i],col=col[i],cex=cex)
	        }
	    else if (mark.time==TRUE && any(deaths==zero.one)) {
		points(xx[deaths==zero.one], 
		       yy[deaths==zero.one],
		       pch=mark[i],col=col[i],cex=cex)
	        }

	    xend <- c(xend,max(xx))
	    yend <- c(yend,min(yy))

	    if (conf.int==TRUE && !is.null(supper)) {
		if (ncurve==1) lty[i] <- lty[i] +1
		yy <- c(firsty, supper[who])
		lines(dostep(xx,yy), lty=lty[i], col=col[i], lwd=lwd[i])
		yy <- c(firsty, slower[who])
		lines(dostep(xx,yy), lty=lty[i], col=col[i], lwd=lwd[i])
	        }
	    }
        }



    invisible(list(x=xend, y=yend))
    }







