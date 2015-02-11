#
# This function is simply an alias for "survfit".  In the Cox model
#  case users often look for the words "baseline hazard"
#
basehaz <- function (fit, centered = TRUE) 
{
    if (!inherits(fit, "coxph")) 
        stop("must be a coxph object")
    sfit <- survfit(fit, se.fit=FALSE)   

    if (!centered) {
        # The right thing to do here is to call survfit with a vector of
        #  all zeros for the "subject to predict".  But if there is a factor
        #  in the model, there may be no subject at all who will give all
        #  zeros, so we post process instead
        zcoef <- ifelse(is.na(coef(fit)), 0, coef(fit))
        offset <- sum(fit$means * zcoef)
        chaz <- sfit$cumhaz * exp(-offset)
    }
    else chaz <- sfit$cumhaz

    new <- data.frame(hazard=chaz, time=sfit$time)

    strata <- sfit$strata
    if (!is.null(strata)) 
        new$strata <- factor(rep(names(strata), strata), levels = names(strata))
    new
}
