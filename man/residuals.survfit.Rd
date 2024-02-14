\name{residuals.survfit}
\alias{residuals.survfit}
\title{IJ residuals from a survfit object.}
\description{
  Return infinitesimal jackknife residuals from a survfit object, for
  the survival, cumulative hazard, or restricted mean time in state (RMTS).
}
\usage{
\method{residuals}{survfit}(object, times, 
    type="pstate", collapse, weighted=FALSE,
    method=1, ...)
}
\arguments{
  \item{object}{a \code{survfit} object}
  \item{times}{a vector of times at which the residuals are desired}
  \item{type}{the type of residual, see below}
  \item{collapse}{add the residuals for all subjects in a cluster. If
    the \code{survfit} object used an \code{id} statement, the default is
    to collapse over that variable.}
  \item{weighted}{weight the residuals by each observation's weight}
  \item{method}{controls a choice of algorithm.  Current an internal
    debugging option.}
  \item{...}{arguments for other methods}
}
\details{
  This function is designed to efficiently compute the residuals at a
small number of time points, also known as the infinitesmal jackknife
(IJ).  Primary usages are the creation of
pseudo-values and IJ estimates of variance.
If the residuals at all time points are needed, e.g. to compute a robust
pointwise confidence interval for the survival curve, then this can be done
more efficiently
using the \code{influence} argument of the underlying \code{survfit}
function.  But be aware that such matrices can get very large. 

The residuals are the impact of each observation or cluster on the
resulting probability in state curves at the given time points,
the cumulative hazard curv at those time points,
or the expected sojourn time in each state up to the given time points.
For a simple Kaplan-Meier the \code{survfit} object contains only the
probability in the "initial" state, i.e., the survival fraction.
For the KM case the sojourn time, the expected amount of time spent in
the initial state, up to the specified endpoint, is more commonly known as the
restricted mean survival time (RMST).
For a multistate model this same quantity is also referred to as the
restricted mean time in state (RMTS).
It can be computed as the area under the respective probability in state curve.
The program allows any of \code{pstate}, \code{surv}, \code{cumhaz},
\code{chaz}, \code{sojourn}, \code{rmst}, \code{rmts} or \code{auc}
for the type argument, ignoring upper/lowercase, so
users can choose whichever abbreviation they like best.

When \code{collapse=TRUE} the result has the cluster identifier (which
defaults to the \code{id} variable) as the dimname for the first dimension.
If the \code{fit} object contains more than one curve, and the same
identifier is reused in two different curves this approach does not work
and the routine will stop with an error.
In principle this is not necessary, e.g., the result could contain two rows
with the same label, showing the separate effect on each curve,
but this was deemed too confusing. 
}

\value{A matrix or array with one row per observation or cluster, and one column
  for each value in \code{times}.  For a multi-state model the three
  dimensions are observation, time and state.  For cumulative hazard,
  the last dimension is the set of transitions.  (A competing risks
  model for instance has 3 states and 2 transitions.)
}
\seealso{\code{\link{survfit}}, \code{\link{survfit.formula}}
}
\examples{
fit <- survfit(Surv(time, status) ~ x, aml)
resid(fit, times=c(24, 48), type="RMTS")
}
% \keyword{ survival }

