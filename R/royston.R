# The D measure proposed in Roysten and Sauerbrei, Stat in Med, 2004
# If eta is the risk score from a Cox model, and v = var(eta), this is
#  essentially R = v/(v+1)
# But not quite: they re-fit the model using a normalized version of r.
#
royston <- function(fit, newdata, ties=TRUE, adjust=FALSE) {
    pi <- 3.141592653589793    # a user might have overwritten the constant 
    if (!inherits(fit, "coxph")) stop("function defined only for coxph models")
    if (inherits(fit, "coxphms")) stop("not defined for multi-state models")

    if (missing(newdata)) {
        eta <- predict(fit)
        y2 <- fit$y
    }
    else {
        eta <- predict(fit, newdata)
        yform <- formula(fit)
        yform[3] <- 1   
        y2 <- model.response(model.frame(yform, data=newdata))
    }
    n <- length(eta)

    R.pm <- var(eta)/(pi^2/6 +var(eta))  # the measure of Kent and O'Quigley
    # eta = X beta, the linear predictor
    # They replace it with a "nicer" one, z = normal scores
    #  If there are ties in eta, replace each with the average normal
    # score to which it matches
    if (ties && any(duplicated(eta))) {
        z <- qnorm((1:n - 3/8)/(n + .25))
        
        # per the paper, take the average z over any ties
        index1 <- match(eta, sort(unique(eta)))
        index2 <- rank(eta, ties.method='first')
        z2 <- tapply(z[index2], index1, mean)
        qhat <- z2[index1]
    }
    else qhat <- qnorm((rank(eta) -3/8)/(n + .25)) # simple case of no ties

    # Do a Cox regression with this "nice" predictor, which has mean 0 and std 1
    #  beta from the Cox model = std of the linear predictor
    rfit <- coxph(y2 ~ qhat)
    beta <- unname(coef(rfit))
    D    <- beta * sqrt(8/pi)  # They consider D the main result
    se.D <- sqrt(rfit$var[1,1]* 8/pi)
    R2   <- beta^2/ (pi^2/6 + beta^2)  # most users will look at R-squared
    R.I  <- beta^2/ (1 + beta^2)       # their R_I statistic

    if (adjust) {
        # overfitting adjustment, related to events per coefficient
        #  if user follows the rule of 10-20 per coef it will often be small
        r <- fit$nevent/(fit$nevent- length(coef(fit)))  # will be > 1
        temp <- (1 + beta^2 -r) /r  # negative values are rare
        D <- sign(beta)*sign(temp) *sqrt(abs(temp) *8/pi)
        se.D <- se.D * abs(beta)/(r*sqrt(abs(temp)))
        R2 <- 1- r*(1-R.I)
    }       
    c(D  = D, "se(D)" = se.D, R.D = R2, R.PM=R.pm)  # return vector
}
