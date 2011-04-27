#
# Reconstruct the Cox model data.  This is done in so many routines
#  that I extracted it out.
# Newer routines use model.matrix.coxph and model.frame.coxph methods.
#
# The "stratax" name is to avoid conflicts with the strata() function, but
#   still allow users to type "strata" as an arg.
#
coxph.getdata <- function(fit, y=TRUE, x=TRUE, stratax=TRUE, offset=FALSE) {
    ty <- fit[['y']]  #avoid grabbing this by accident due to partial matching
    tx <- fit[['x']]  #  for x, fit$x will get fit$xlevels --> not good
    strat <- fit$strata
    Terms <- fit$terms
    if (is.null(attr(Terms, 'offset'))) offset <- FALSE
    if (offset) x<- TRUE
    if (!inherits(Terms, 'terms'))
	    stop("invalid terms component of fit")
    strats <- attr(Terms, "specials")$strata
    if (length(strats)==0) stratax <- FALSE

    if ( (y && is.null(ty)) || (x && is.null(tx)) ||
	     (stratax && is.null(strat)) || offset) {
	# get the model frame
	m <- model.frame(fit)

	# Pull things out
	if (y && is.null(ty)) ty <- model.extract(m, 'response')

	if (offset) toff <- model.extract(m, 'offset')

	# strata was saved in the fit if and only if x was
	if ((x || stratax) && is.null(tx)) {
	    if (stratax) {
		temp <- untangle.specials(Terms, 'strata', 1)
		strat <- strata(m[temp$vars], shortlabel=T)
		}
	    if (x) tx <- model.matrix(fit, data=m)
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
