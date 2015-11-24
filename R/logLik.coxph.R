#
# The AIC function depends on a logLik method
#
logLik.coxph <- function(object, ...) {
    out <- object$loglik[2]
    if (!is.null(object$df)) attr(out, "df") <- object$df[2]
    else  attr(out, 'df') <- sum(!is.na(coefficients(object)))
    attr(out, "nobs") <- object$nevent
    class(out) <- 'logLik'
    out
    }

# Cox models with no covariates
logLik.coxph.null <- function(object, ...) {
    out <- object$loglik[1]
    attr(out, "df") <- 0
    attr(out, "nobs") <- object$nevent
    class(out) <- "loglik"
    out
}
    
logLik.survreg <- function(object, ...) {
    out <- object$loglik[2]
    dd <- diag(object$var)
    if (!is.null(object$df)) attr(out, "df") <- sum(object$df)
    else attr(out, 'df') <- sum(!is.na(dd) & dd > 0)
#    attr(out, "nobs") <- sum(object$df) + object$df.residual
    class(out) <- 'logLik'
    out
    }
