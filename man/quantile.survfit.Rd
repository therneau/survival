\name{quantile.survfit}
\alias{quantile.survfit}
\title{Quantiles from a survfit object}
\description{Retrieve quantiles and confidence intervals for them from
  a survfit object.
}
\usage{
\method{quantile}{survfit}(x, probs = c(0.25, 0.5, 0.75), conf.int = TRUE,
...)
}
\arguments{
  \item{x}{a result of the survfit function}
  \item{probs}{numeric vector of probabilities with values in [0,1].}
  \item{conf.int}{should lower and upper confidence limits be returned?}
  \item{...}{optional arguments for other methods}
}
\details{
  The kth quantile for a survival curve S(t) is the location at which
  a horizontal line at height 1-k intersects the plot of S(t).
  Since S(t) is a step function, it is possible for the curve to have a
  horizontal segment at exactly 1-k, in which case the midpoint of the
  horizontal segment is returned.  This mirrors the standard behavior of
  the median when data is uncensored.  If the survival curve does not
  fall to 1-k, then that quantile is undefined.

  Confidence limits for the values are based on the intersection of the
  horizontal line at 1-k with the upper and lower limits for the
  survival curve.  Hence confidence limits use the same
  p-value as was in effect when the curve was created, and will differ
  depending on the \code{conf.type} option of \code{survfit}.
  If the survival curves have no confidence bands, confidence limits for
  the quantiles are not available.
}
\value{
    The quantiles will be a vector if the \code{survfit} object contains
    only a single curve, otherwise it will be a matrix or array.  In
    this latter case the second dimension will index the quantiles.
    
    If confidence limits are requested, then result will list with components
    \code{quantile}, \code{lower}, and \code{upper}, otherwise the
    vector or matrix of quantiles.
}

\author{Terry Therneau}
\seealso{\code{\link{survfit}}, \code{\link{print.survfit}},
  \code{\link{qsurvreg}}
}
\examples{
fit <- survfit(Surv(time, status) ~ ph.ecog, data=lung)
quantile(fit)
}
\keyword{ survival }
