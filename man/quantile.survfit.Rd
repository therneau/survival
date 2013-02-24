\name{quantile.survfit}
\alias{quantile.survfit}
\alias{quantile.survfitms}
\title{Quantiles from a survfit object}
\description{Retrieve quantiles and confidence intervals for them from
  a survfit object.
}
\usage{
\method{quantile}{survfit}(x, probs = c(0.25, 0.5, 0.75), conf.int = TRUE,
  tolerance= sqrt(.Machine$double.eps), ...)
\method{quantile}{survfitms}(x, probs = c(0.25, 0.5, 0.75), conf.int = TRUE,
  tolerance= sqrt(.Machine$double.eps), ...)
}
\arguments{
  \item{x}{a result of the survfit function}
  \item{probs}{numeric vector of probabilities with values in [0,1]}
  \item{conf.int}{should lower and upper confidence limits be returned?}
  \item{tolerance}{tolerance for checking that the suvival curve exactly
    equals one of the quantiles}
  \item{...}{optional arguments for other methods}
}
\details{
  The kth quantile for a survival curve S(t) is the location at which
  a horizontal line at height p= 1-k intersects the plot of S(t).
  Since S(t) is a step function, it is possible for the curve to have a
  horizontal segment at exactly 1-k, in which case the midpoint of the
  horizontal segment is returned.  This mirrors the standard behavior of
  the median when data is uncensored.  If the survival curve does not
  fall to 1-k, then that quantile is undefined.

  In order to be consistent with other quantile functions, the argument
  \code{prob} of this function applies to the cumulative distribution
  function F(t) = 1-S(t).  

  Confidence limits for the values are based on the intersection of the
  horizontal line at 1-k with the upper and lower limits for the
  survival curve.  Hence confidence limits use the same
  p-value as was in effect when the curve was created, and will differ
  depending on the \code{conf.type} option of \code{survfit}.
  If the survival curves have no confidence bands, confidence limits for
  the quantiles are not available.

  When a horizontal seqment of the survival curve exactly matches one of
  the requested quantiles the returned value will be the midpoint of the
  horizontal seqment; this agrees with the usual definition of a median
  for uncensored data.  Since the survival curve is computed as a series
  of products, however, there may be round off error.
  Assume for instance a sample of size 20 with no tied times and no
  censoring.  The survival curve after the 10th death is
  (19/20)(18/19)(17/18) ... (10/11) = 10/20, but the computed result will
  not be exactly 0.5. Any horizontal segment whose absolute difference
  with a requested percentile is less than \code{tolerance} is
  considered to be an exact match.
}
\value{
    The quantiles will be a vector if the \code{survfit} object contains
    only a single curve, otherwise it will be a matrix or array.  In
    this case the last dimension will index the quantiles.
    
    If confidence limits are requested, then result will be a list with
    components
    \code{quantile}, \code{lower}, and \code{upper}, otherwise it is the
    vector or matrix of quantiles.
}

\author{Terry Therneau}
\seealso{\code{\link{survfit}}, \code{\link{print.survfit}},
  \code{\link{qsurvreg}}
}
\examples{
fit <- survfit(Surv(time, status) ~ ph.ecog, data=lung)
quantile(fit)

cfit <- coxph(Surv(time, status) ~ age + strata(ph.ecog), data=lung)
csurv<- survfit(cfit, newdata=data.frame(age=c(40, 60, 80)),
                  conf.type ="none")
temp <- quantile(csurv, 1:5/10)
temp[2,3,]  # quantiles for second level of ph.ecog, age=80
quantile(csurv[2,3], 1:5/10)  # quantiles of a single curve, same result
}
\keyword{ survival }
