\name{residuals.survfit}
\alias{residuals.survfit}
\title{IJ residuals from a survfit object.}
\description{
  Obtain infinitesimal jackknife residuals from a survfit object, for
  the survival, cumulative hazard, or restricted mean time in state (RMTS).
}
\usage{
\method{residuals}{survfit}(object, times, type = c("surv", "RMTS", "cumhaz"),
  weighted=FALSE, collapse=FALSE)
}
\arguments{
  \item{object}{a \code{survfit} object}
  \item{times}{a vector of times at which the residuals are desired}
  \item{type}{the type of residual}
  \item{weighted}{weight the residuals by each observation's weight}
  \item{collapse}{add the residuals for all subjects in a cluster}
}
\details{
This function is designed to efficiently compute the residuals at a
small number of time points; a primary use is the creation of
pseudo-values.
If the residuals at all time points are needed, e.g. to compute a robust
pointwise confidence interval for the survival curve, this can be done
using the \code{influence} argument of the underlying \code{survfit}
function.  (But be aware that such matrices can get very large).

For an ordinary alive->death survival curve, mean time in state is only
of interest for the "alive" state, and in this case the acronym
"restricted mean survival time" (RMST) is much more common than RMTS,
which is more general.  
}
\value{A matrix or array with one row per observation or cluster, and one column
  for each value in \code{times}.  In a multi-state model the state will
  be the third dimension, e.g., restricted mean time in each state.
}
\seealso{\code{\link{survfit}}, \code{\link{survfit.formula}}
}
\examples{
fit <- survfit(Surv(time, status) ~ x, aml)
resid(fit, times=c(24, 48), type="RMTS")
}
% \keyword{ survival }

