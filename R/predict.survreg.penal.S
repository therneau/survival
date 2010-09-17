# $Id: predict.survreg.penal.S 11166 2008-11-24 22:10:34Z therneau $
#
# This routine just stops disastrous arithmetic for models with sparse
# terms.  A placeholder until the proper sparse terms actions are inserted.
#
predict.survreg.penal <- function(object, ...) {
    pterms <- object$pterms
    if (any(pterms==2))
	    stop("Predictions not available for sparse models")
    NextMethod('predict')
    }
