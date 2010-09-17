# $Id: lines.survfit.coxph.S 11166 2008-11-24 22:10:34Z therneau $
lines.survfit.coxph <- function(x, mark.time=FALSE, ...) {
    if (is.logical(mark.time) & mark.time)
	    stop("Invalid value for mark.time")
    invisible(NextMethod('lines', mark.time=mark.time))
    }
