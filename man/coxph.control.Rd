\name{coxph.control}
\alias{coxph.control}
\title{Ancillary arguments for controling coxph fits}
\description{
This is used to set various numeric parameters controling a Cox model fit.
Typically it would only be used in a call to \code{coxph}.
}
\usage{
coxph.control(eps = 1e-09, toler.chol = .Machine$double.eps^0.75,
iter.max = 20, toler.inf = sqrt(eps), outer.max = 10)
}
\arguments{
  \item{eps}{Iteration continues until the relative change in the log partial
likelihood is less than eps.  Must be positive.}
  \item{toler.chol}{Tolerance for detection of singularity during a Cholesky
decomposion of the variance matrix, i.e., for detecting a redundant predictor
variable.}
  \item{iter.max}{Maximum number of iterations to attempt for convergence.}
  \item{toler.inf}{Tolerance criteria for the warning message about a possible
infinite coefficient value.}
  \item{outer.max}{For a penalized coxph model, e.g. with pspline terms, there
is an outer loop of iteration to determine the penalty parameters; maximum
number of iterations for this outer loop.}
}

\value{
  a list containing the values of each of the above constants
}
\author{Terry Therneau }
\seealso{\code{\link{coxph}}
  }
\keyword{survival}
