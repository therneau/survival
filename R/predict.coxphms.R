predict.coxphms <- function(object, newdata, 
                       type=c("lp", "risk", "expected", "terms", "survival"),
                       se.fit=FALSE, na.action=na.pass,
                       terms=names(object$assign), collapse, 
                       reference=c("strata", "sample"), ...) {
    if (missing(newdata) && type %in% c("lp", "risk")) NextMethod();
    else stop("predict method not yet available for multistate coxph")
}

residuals.coxphms <- function(object, type=c("martingale", "deviance", "score",
                                             "schoenfeld",
			  "dfbeta", "dfbetas", "scaledsch","partial"),
                          collapse=FALSE, weighted=FALSE, ...) {
    stop("residuals method not yet available for multistate coxph")
}


