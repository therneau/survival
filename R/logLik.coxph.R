#
# The AIC function depends on a logLik method
#
logLik.coxph <- function(object, ...) {
    out <- object$loglik[2]
    attr(out, 'df') <- sum(!is.na(coefficients(object)))
    class(out) <- 'logLik'
    out
    }

logLik.survreg <- function(object, ...) {
    out <- object$loglik[2]
    dd <- diag(object$var)
    attr(out, 'df') <- sum(!is.na(dd) & dd > 0)
    attr(out, "nobs") <- object$df + object$df.residual
    class(out) <- 'logLik'
    out
    }
