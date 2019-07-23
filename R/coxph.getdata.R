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

    if (is.null(fit$call$weights)) twt <- rep(1.0, fit2$n)
    else twt<- fit[["weights"]]

    strat <- fit$strata
    Terms <- fit$terms
    if (is.null(attr(Terms, 'offset'))) offset <- FALSE
    if (offset) x<- TRUE  # can't grab offset without x
    if (!inherits(Terms, 'terms'))
	    stop("invalid terms component of fit")
    strats <- attr(Terms, "specials")$strata
    if (length(strats)==0) stratax <- FALSE

    if ( (y && is.null(ty)) || (x && is.null(tx)) || 
         (weights && is.null(twt)) ||
	     (stratax && is.null(strat)) || offset || weights) {
	# get the model frame
	mf <- stats::model.frame(fit)

	# Pull things out
	if (y && is.null(ty)) ty <- model.extract(mf, "response")
        if (is.null(twt))     twt<- model.extract(mf, "weights")
	if (offset) toff <- model.extract(mf, 'offset')

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
            if (offset) toff <- offset[xstack$rindex]

            # And last, toss missing values, which had been deferred
            ismiss <- is.nan(ty) | apply(is.na(tx), 1, any)
            if (offset) {
                ismiss <- ismiss | is.nan(offset)
                offset <- offset[imiss]
            }       
	    if (any(ismiss)) {
                if (y) ty<- ty[!ismiss]
                tx <- tx[!ismiss,,drop=FALSE]
                strat <- strat[!ismiss]
            }       
	}
    }
    else if (offset)
       toff <- fit$linear.predictors -(c(tx %*% fit$coef) - 
                                        sum(fit$means*fit$coef))

    temp <- list()
    if (y) temp$y <- ty
    if (x) temp$x <- tx
    if (stratax)  temp$strata <- strat
    if (offset)  temp$offset <- toff
    temp
    }
