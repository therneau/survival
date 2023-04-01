\name{brier}
\alias{brier}
\title{Compute the Brier score for a Cox model}
\description{Compute the Brier score, for a coxph model}
\usage{
brier(fit, times, newdata, ties = TRUE, detail = FALSE, timefix = TRUE, 
      efron = FALSE)
}
\arguments{
  \item{fit}{result of a \code{coxph} fit}
  \item{times}{time points at which to create the score}
  \item{newdata}{optional, used to validate a prior fit with new data}
  \item{ties}{if TRUE, treate tied event/censoring times properly}
  \item{detail}{if TRUE, the returned object has more detail.  This can
  be useful for debugging or for instruction.}
  \item{timefix}{deal with near ties in the data.  See the tied times vignette.}
  \item{efron}{use the same survival estimate for the NULL model as was
    used in the coxph call}
}
\details{
Far more details are found in the vignette. At any time point tau, the
scaled Brier score is essentially the R-squared statistic where y =
the 0/1 variable "event at or before tau", yhat is the probability of an
event by tau, as predicted by the model, and the ybar is the predicted
probablity without covariate, normally from a Kaplan-Meier.
If \eqn{R^2= 1- \sum(y- \hat y)^2/\sum (y- \mu)^2}{%
R-squared= 1- mse(fitted)/mse(null)}, the Brier score is formally only
the numerator of the second term.  The rescaled value is much more
useful, however.

Many, perhaps even most algorithms do not properly deal with a tied
censoring time/event time pair.  The \code{tied} option is present
mostly verify that we get the same answer, when we make the same
mistake.  The numerical size of the inaccuracy is very small; just large
enough to generate concern that this function is incorrect.

A sensible argument can be made that the NULL model should be a
\code{coxph} call with no covariates, rather than the Kaplan-Meier;
but it turns out that the effect is very slight.
This is allowed by the \code{efron} argument.
}
\value{ a list with components
  \item{rsquared}{the \eqn{R^2}{R-squared} value, a scaled Brier score.
  This will be a vector with one entry for each time point.}
\item{brier}{the brier score, a vector with one entry per time point}
\item{times}{the time points at which the score was computed}
}

\author{Terry Therneau}
\seealso{\code{\link{rttright}}}
\examples{
cfit <- coxph(Surv(rtime, recur) ~ age + meno + size + pmin(nodes,11), 
              data= rotterdam)
round(cfit$concordance["concordance"], 3)  # some predictive power
brier(cfit, times=c(4,6)*365.25)   # values at 4 and 6 years
}
\keyword{ survival }
