#
# This function is simply an alias for "survfit".  In the Cox model
#  case users often look for the words "baseline hazard"
#
basehaz <- function (fit, centered = TRUE) 
{
    if (!inherits(fit, "coxph")) 
        stop("must be a coxph object")
    if (!centered) {
        tdata <- as.list(0*fit$means)  #dummy data set of all zeros
        sfit <- survfit(fit, newdata=tdata, se.fit=FALSE)
    }
    else sfit <- survfit(fit, se.fit=FALSE)

    new <- data.frame(hazard=sfit$cumhaz, time=sfit$time)

    strata <- sfit$strata
    if (!is.null(strata)) 
        new$strata <- factor(rep(names(strata), strata), levels = names(strata))
    new
}
