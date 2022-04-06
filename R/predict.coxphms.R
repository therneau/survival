predict.coxphms <- function(object, newdata, 
                       type=c("lp", "risk", "expected", "terms", "survival"),
                       se.fit=FALSE, na.action=na.pass,
                       terms=names(object$assign), collapse, 
                       reference=c("strata", "sample"), ...) {

    type <- match.arg(type)
    if (missing(newdata) && (type %in% c("lp", "risk"))) NextMethod()
    else stop("predict method not yet available for multistate coxph")
}


