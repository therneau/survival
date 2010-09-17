#
# The AIC function depends on a logLik method
#
logLik.coxph <- function(x) {
    out <- x$loglik[2]
    attr(out, 'df') <- sum(!is.na(coefficients(x)))
    class(out) <- 'logLik'
    out
    }

logLik.survreg <- function(x) {
    out <- x$loglik[2]
    attr(out, 'df') <- sum(diag(fit$var) > 0)
    class(out) <- 'logLik'
    out
    }
