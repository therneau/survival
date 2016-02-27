# A primary purpose of this function is to not print out everything,
#  e.g. the terms structure

print.pyears <- function(x, ...) {
    if (!is.null(cl<- x$call)) {
	cat("Call: ")
	dput(cl)
	cat("\n")
	}
    cat("n=", x$observations)
    if (length(x$na.action))
	cat(" (", naprint(x$na.action), ")\n", sep="")
    else cat("\n")
    if (length(x$offtable) && x$offtable>0)
        cat("  ", format(x$offtable), " years off table\n\n")
    else cat("\n")

    ii <- match(c("event", "pyears", "expected", "n"), names(x), nomatch=0)
    print(x[ii])
    invisible(x)
}
