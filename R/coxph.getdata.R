#
# Reconstruct the Cox model data.  This is done in many routines.
# Users use model.matrix.coxph and model.frame.coxph methods, but they
#  do not extract strata or offset.
#
# The "stratax" name is to avoid conflicts with the strata() function, but
#   still allow users to type "strata" as an arg.
#
coxph.getdata <- function(fit, y=TRUE, x=TRUE, stratax=TRUE, 
                          weights=TRUE, offset=FALSE) {
    ty <- fit[['y']]  #avoid grabbing this by accident due to partial matching
    tx <- fit[['x']]  #  for x, fit$x will get fit$xlevels --> not good

    # if x or y is present, use it to set n
    if (!is.null(ty)) n <- nrow(ty)
    else if (!is.null(tx)) n <- nrow(tx)
    else n <- NULL

    Terms <- fit$terms
    if (!inherits(Terms, 'terms'))
	    stop("invalid terms component of fit")

    # Avoid calling model.frame unless we have to: fill in weights and/or
    #  offset when they were not present
    if (!is.null(n)) {
        if (is.null(fit$call$weights)) twt <- rep(1,n)
        if (is.null(attr(terms(fit), "offset"))) toff <- rep(0, n)
    }
    else {
        toff <- NULL
        twt<- fit[["weights"]]
    }

    strats <- attr(Terms, "specials")$strata
    if (length(strats)==0) stratax <- FALSE
    strat <- fit$strata

    if ((y && is.null(ty)) || (x && is.null(tx)) || 
        (weights && is.null(twt)) ||
	(stratax && is.null(strat)) || (offset && is.null(toff))) {
	# get the model frame
	mf <- stats::model.frame(fit)
        n <- nrow(mf)
        
	# Pull things out
        if (is.null(twt) && weight) {
            twt <- model.extract(mf, "weights")
            if (is.null(twt)) twt <- rep(1.0, n)
        }
            
	if (offset && is.null(toff)) {
            toff <- model.extract(mf, 'offset')
            if (is.null(toff)) toff <- rep(0.0, n)
        }

	if (y && is.null(ty)) ty <- model.extract(mf, "response")

	# strata was saved in the fit if and only if x was
	if ((x || stratax) && is.null(tx)) {
	    if (stratax) {
		temp <- untangle.specials(Terms, 'strata', 1)
		strat <- strata(mf[temp$vars], shortlabel=T)
		}
	    tx <- model.matrix.coxph(fit, data=mf)
            }

        if (inherits(fit, "coxphms")) {
            # Expand the x matrix.  First recreate istate
            id <- model.extract(mf, "id")
            istate <- model.extract(mf, "istate")
            check <- survcheck2(ty, id, istate)
            
            # Now expand the data
            xstack <- stacker(fit$cmap, check$istate, tx, ty, 
                              as.integer(strata), check$states)
            tx <- xstack$X
            ty <- xstack$Y
            strat <- xstack$strata
            stratax <- TRUE
            if (offset) toff <- toff[xstack$rindex]
            if (weight) twt  <- twt[xstack$rindex]

            # And last, toss missing values, which had been deferred
            ismiss <- is.nan(ty) | apply(is.na(tx), 1, any)
            if (offset) ismiss <- ismiss | is.nan(toff)
            if (weight) ismiss <- ismiss | is.nan(twt)

            if (any(ismiss)) {
                if (offset) toff <- toff[!ismiss]
                if (weight) twt  <- twt[!ismiss]
            
                if (y) ty<- ty[!ismiss]
                if (x) tx <- tx[!ismiss,,drop=FALSE]
                if (stratax) strat <- strat[!ismiss]
            }       
	}
    }

    temp <- list()
    if (y) temp$y <- ty
    if (x) temp$x <- tx
    if (stratax)  temp$strata <- strat
    if (offset)  temp$offset <- toff
    if (weights) temp$weights <- twt
    temp
    }
