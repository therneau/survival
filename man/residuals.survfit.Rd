\name{residuals.survfit}
\alias{residuals.survfit}
\title{IJ residuals from a survfit object.}
\description{
  Return infinitesimal jackknife residuals from a survfit object, for
  the survival, cumulative hazard, or restricted mean time in state (RMTS).
}
\usage{
\method{residuals}{survfit}(object, times, 
    type, collapse=TRUE, weighted=FALSE, method=1, ...)
}
\arguments{
  \item{object}{a \code{survfit} object}
  \item{times}{a vector of times at which the residuals are desired}
  \item{type}{the type of residual, see below}
  \item{collapse}{add the residuals for all subjects in a cluster}
  \item{weighted}{weight the residuals by each observation's weight}
  \item{method}{controls a choice of algorithm.  Current an internal
    debugging option.}
  \item{...}{arguments for other methods}
}
\details{
This function is designed to efficiently compute the residuals at a
small number of time points; a primary use is the creation of
pseudo-values.
If the residuals at all time points are needed, e.g. to compute a robust
pointwise confidence interval for the survival curve, then this can be done
more efficiently
using the \code{influence} argument of the underlying \code{survfit}
function.  But be aware that such matrices can get very large. 

The residuals are the impact of each observation or cluster on the
resulting \code{surv}, \code{pstate} or code{cumhaz} component of a
\code{survfit} object.  The restricted mean survival time (RMST)
is E[min(survival, t)], and can be computed as the area under the
survival curve.
A multi-state survival curve has a \code{pstate} component instead of
\code{surv}, and the area under each component is the restricted mean time
in state (RMTS).
For simplicity, in this function \code{surv} and \code{pstate} are
treated as synonyms for the \code{type} argument,
and likewise for \code{rmst}, \code{rmts} or their capitalized versions.

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
  dimensions are observation, state, and time.
}
\seealso{\code{\link{survfit}}, \code{\link{survfit.formula}}
}
\examples{
fit <- survfit(Surv(time, status) ~ x, aml)
resid(fit, times=c(24, 48), type="RMTS")
}
% \keyword{ survival }

