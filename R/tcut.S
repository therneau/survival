# $Id: tcut.S 11166 2008-11-24 22:10:34Z therneau $
tcut <-  function (x, breaks, labels, scale=1){
    # avoid some problems with dates
    x <- as.numeric(x)
    breaks <- as.numeric(breaks)

    if(length(breaks) == 1) {
	if(breaks < 1)
		stop("Must specify at least one interval")
	if(missing(labels))
		labels <- paste("Range", seq(length = breaks))
	else if(length(labels) != breaks)
		stop("Number of labels must equal number of intervals")
	r <- range(x[!is.na(x)])
	r[is.na(r)] <- 1
	if((d <- diff(r)) == 0) {
		r[2] <- r[1] + 1
		d <- 1
	    }
	breaks <- seq(r[1] - 0.01 * d, r[2] + 0.01 * d, length = breaks +1)
	}
    else {
	if(is.na(adb <- all(diff(breaks) >= 0)) || !adb)
	   stop("breaks must be given in ascending order and contain no NA's")
	if(missing(labels))
	    labels <- paste(format(breaks[ - length(breaks)]),
			"+ thru ", format(breaks[-1]), sep = "")
	else if(length(labels) != length(breaks) - 1)
	   stop("Number of labels must be 1 less than number of break points")
	}

    temp <- structure(x*scale, cutpoints=breaks*scale, labels=labels)
    if (is.R()) class(temp) <- 'tcut'
    else        oldClass(temp) <- 'tcut'
    temp
    }

"[.tcut" <- function(x, ..., drop=FALSE) {
    atts <- attributes(x)
    x <- unclass(x)[..1]
    attributes(x) <- atts
    if (is.R()) class(x) <- 'tcut'
    else        oldClass(x) <- 'tcut'
    x
    }
