#
# new generic methods sprout up regularly
# This is a convenient place to add them
#

# The AER package has a deviance function for tobit models, but I don't think
#  that deviance is a useful quantity there.
#
fitted.survreg <- function(object, ...)
  predict(object, type = "response", se.fit = FALSE)

fitted.coxph <- function(object, ...) object$linear.predictors

# The nobs result is intended for computing BIC; for coxph we think that
#  this is the number of events.
nobs.survreg <- function(object, ...)
  length(object$linear.predictors)
nobs.coxph <- function(object, ...) object$nevent

weights.survreg <- function(object, ...)
  model.weights(model.frame(object))

weights.coxph <- function(object, ...) {
    if (!is.null(object$weights)) object$weights
    else model.weights(model.frame(object))
}
