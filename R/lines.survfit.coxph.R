lines.survfit.coxph <- function(x, mark.time=FALSE, ...) {
    if (is.logical(mark.time) & mark.time)
	    stop("Invalid value for mark.time")
    invisible(NextMethod('lines', mark.time=mark.time))
    }
