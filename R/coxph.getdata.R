#
# Reconstruct the Cox model data.  This is done in many routines.
# Users use model.matrix.coxph and model.frame.coxph methods, but they
#  do not extract strata, offset, or istate
#
# The "stratax" name is to avoid conflicts with the strata() function, but
#   still allow users to type "strata" as an arg.
#
coxph.getdata <- function(fit, y=TRUE, x=TRUE, stratax=TRUE, 
                          weights=TRUE, offset=FALSE, id=TRUE, cluster=TRUE) {
    ty <- fit[['y']]  #avoid grabbing this by accident due to partial matching
    tx <- fit[['x']]  #  for x, fit$x will get fit$xlevels --> not good
    twt <- fit[["weights"]]
    toff <- fit[["offset"]]
    if (is.null(fit$call$id)) id <- FALSE  # there is no id to return
    if (is.null(fit$call$cluster)) cluster <- FALSE

    # if x or y is present, use it to set n
    if (!is.null(ty)) n <- nrow(ty)
    else if (!is.null(tx)) n <- nrow(tx)
    else n <- NULL

    coxms <- inherits(fit, "coxphms")
    Terms <- fit$terms
    if (!inherits(Terms, 'terms'))
	    stop("invalid terms component of fit")

    # Avoid calling model.frame unless we have to: fill in weights and/or
    #  offset when they were not present.  But we can only do it successfully
    #  if we know n.
    if (!is.null(n)) {
        if (is.null(fit$call$weights)) twt <- rep(1,n)
        if (is.null(attr(terms(fit), "offset"))) toff <- rep(0, n)
    }

    strat <- fit$strata
    strats <- attr(Terms, "specials")$strata
    if (length(strats)==0 && length(strat)==0 & !coxms) stratax <- FALSE

    if ((y && is.null(ty)) || (x && is.null(tx)) || 
        (weights && is.null(twt)) ||  cluster || id ||
	(stratax && is.null(strat)) || (offset && is.null(toff)) ||
        !is.null(fit$call$istate)) {
	# get the model frame
	mf <- stats::model.frame(fit)
        n <- nrow(mf)

	# Pull things out
        if (weights) {
            twt <- model.extract(mf, "weights")
            if (is.null(twt)) twt <- rep(1.0, n)
        }
            
	if (offset) {
            toff <- model.extract(mf, 'offset')
            if (is.null(toff)) toff <- rep(0.0, n)
        }
        if (id) idx <- model.extract(mf, "id")
        if (cluster) clusterx <- model.extract(mf, "cluster")

        if (inherits(fit, "coxphms")) {
            # we need to call stacker
            idx <- model.extract(mf, "id")
            istate <- model.extract(mf, "istate")
            ty <- model.response(mf)
            if (is.null(fit$timefix) || fit$timefix) ty <- aeqSurv(ty) 
            check <- survcheck2(ty, idx, istate)
            tx <- model.matrix.coxph(fit, data=mf)
            if (length(strats)) {
  		temp <- untangle.specials(Terms, 'strata', 1)
		strat <- as.integer(strata(mf[temp$vars], shortlabel=T))
            }
            else strat <- NULL
            # Now expand the data
            xstack <- stacker(fit$cmap, fit$smap, as.integer(check$istate), tx, ty, 
                              strat, check$states)
            tx <- xstack$X
            ty <- xstack$Y
            strat <- xstack$strata
            stratax <- TRUE
            if (offset) toff <- toff[xstack$rindex]
            if (weights) twt  <- twt[xstack$rindex]
            if (id) idx <- idx[xstack$rindex]
            if (cluster) clusterx <- clusterx[xstack$rindex]

            # And last, toss missing values, which had been deferred
            ismiss <- is.nan(ty) | apply(is.na(tx), 1, any)
            if (offset) ismiss <- ismiss | is.nan(toff)
            if (weights) ismiss <- ismiss | is.nan(twt)
            if (any(ismiss)) {
                if (offset) toff <- toff[!ismiss]
                if (weights) twt  <- twt[!ismiss]
                if (y) ty<- ty[!ismiss]
                if (x) tx <- tx[!ismiss,,drop=FALSE]
                if (stratax) strat <- strat[!ismiss]
                if (id) idx <- idx[!ismiss]
                if (cluster) clusterx <- clusterx[!ismiss]
            }       
        } 
        else { # not multi-state, or everything was there
            if (y && is.null(ty)) {
                ty <- model.extract(mf, "response")
                if (is.null(fit$timefix) || fit$timefix) ty <- aeqSurv(ty)
            }

            # strata was saved in the fit if and only if x was
            if ((x || stratax) && is.null(tx)) {
                if (stratax) {
                    temp <- untangle.specials(Terms, 'strata', 1)
                    strat <- strata(mf[temp$vars], shortlabel=T)
		}
                tx <- model.matrix.coxph(fit, data=mf)
            }
        }   
    }

    temp <- list()
    if (y) temp$y <- ty
    if (x) temp$x <- tx
    if (stratax)  temp$strata <- strat
    if (offset)  temp$offset <- toff
    if (weights) temp$weights <- twt
    if (id) temp$id <- idx
    if (cluster) temp$cluster <- clusterx
    temp
    }
