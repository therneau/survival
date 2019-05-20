# $Id: labels.survreg.S 11166 2008-11-24 22:10:34Z therneau $
labels.survreg <- function(object, ...)
	attr(object$terms, "term.labels")

