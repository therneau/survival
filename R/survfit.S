survfit <- function (formula, ...) {
    Call <- match.call()
    # Real tricky -- find out if the first arg is "Surv(...)" without
    #  evaluating it. 
    # At one time, we allowed this to act like a formula, by adding the ~1
    #    on for the user.  Very non-standard, and a bad idea.
    # We removed all refernce to this useage from the documentation several
    #    years ago.  The line below is a temporary, "nice" failure message
    #    for this no-longer supported use.  In a couple more years, axe this
    #    section.
    #
    if ((mode(Call[[2]]) == 'call' &&  Call[[2]][[1]] == as.name('Surv'))
		|| inherits(formula, 'Surv'))  {
	stop(paste("Survfit requires a formula or a coxph fit as the",
		   "first argument"))
	}

    UseMethod('survfit', formula)
    }

# The subscript function is bundled in here, although used most
#  often in plotting

"[.survfit" <- function(x, ..., drop=FALSE) {
    if (missing(..1)) i<- NULL  else i <- sort(..1)
    if (missing(..2)) j<- NULL  else j <- ..2
    if (is.null(x$strata)) {
        if (is.matrix(x$surv)) {
            x$surv <- x$surv[,i,drop=drop]
            if (!is.null(x$std.err)) x$std.err <- x$std.err[,i,drop=drop]
            if (!is.null(x$upper)) x$upper <- x$upper[,i,drop=drop]
            if (!is.null(x$lower)) x$lower <- x$lower[,i,drop=drop]
            }
        else warning("Survfit object has only a single survival curve")
        }
    else {
        if (is.null(i)) keep <- seq(along.with=x$time)
        else {
            if (is.character(i)) {
                strat <- rep(names(x$strata), x$strata)
                indx  <- match(i, names(x$strata))
                if (any(is.na(indx))) {
                    stop(paste("subscript(s)", 
                               paste(i[is.na(indx)], collapse=' '),
                               'not matched'))
                    }
                }
            else {
                strat <- rep(1:length(x$strata), x$strata)
                indx <- i
                }
            keep <- seq(along.with=strat)[match(strat, i, nomatch=0)>0]
            if (length(i) <=1) x$strata <- NULL
            else               x$strata  <- x$strata[indx]

            x$n       <- x$n[indx]
            x$time    <- x$time[keep]
            x$n.risk  <- x$n.risk[keep]
            x$n.event <- x$n.event[keep]
            x$n.censor<- x$n.censor[keep]
            if (!is.null(x$enter)) x$enter <- x$enter[keep]
            }
        if (is.matrix(x$surv)) {
            if (is.null(j)) {
                x$surv <- x$surv[keep,,drop=drop]
                if (!is.null(x$std.err)) 
                        x$std.err <- x$std.err[keep,,drop=drop]
                if (!is.null(x$upper)) x$upper <-x$upper[keep,,drop=drop]
                if (!is.null(x$lower)) x$lower <-x$lower[keep,,drop=drop]
                }
            else {
                x$surv <- x$surv[keep,j]
                if (!is.null(x$std.err)) x$std.err <- x$std.err[keep,j]
                if (!is.null(x$upper)) x$upper <- x$upper[keep,j]
                if (!is.null(x$lower)) x$lower <- x$lower[keep,j]
                }
            }
        else {
            x$surv <- x$surv[keep]
            if (!is.null(x$std.err)) x$std.err <- x$std.err[keep]
            if (!is.null(x$upper)) x$upper <- x$upper[keep]
            if (!is.null(x$lower)) x$lower <- x$lower[keep]
            }
        }
    x
    }
