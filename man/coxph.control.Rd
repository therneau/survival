\name{coxph.control}
\alias{coxph.control}
\title{Ancillary arguments for controlling coxph fits}
\description{
This is used to set various numeric parameters controlling a Cox model fit.
Typically it would only be used in a call to \code{coxph}.
}
\usage{
coxph.control(eps = 1e-09, toler.chol = .Machine$double.eps^0.75,
iter.max = 20, toler.inf = sqrt(eps), outer.max = 10, timefix=TRUE)
}
\arguments{
  \item{eps}{Iteration continues until the relative change in the log partial
    likelihood is less than eps, or the absolute change is less than
    sqrt(eps).  Must be positive.}
  \item{toler.chol}{Tolerance for detection of singularity during a Cholesky
decomposition of the variance matrix, i.e., for detecting a redundant predictor
variable.}
  \item{iter.max}{Maximum number of iterations to attempt for convergence.}
  \item{toler.inf}{Tolerance criteria for the warning message about a possible
infinite coefficient value.}
  \item{outer.max}{For a penalized coxph model, e.g. with pspline terms, there
is an outer loop of iteration to determine the penalty parameters; maximum
number of iterations for this outer loop.}
\item{timefix}{Resolve any near ties in the time variables.}
}

\value{
  a list containing the values of each of the above constants
}
\details{
  The convergence tolerances are a balance.  Users think they want THE maximum
  point of the likelihood surface, and for well behaved data sets where
  this is quadratic near the max a high accuracy is fairly inexpensive:
  the number of correct digits approximately doubles with each iteration.
  Conversely, a drop of .0001 from the maximum in any given direction
  will be correspond to only about 1/20 of a standard error change in the
  coefficient.  Statistically, more precision than this is straining at a
  gnat.  Based on this the author originally had set the tolerance to
  1e-5, but relented in the face of multiple
  "why is the answer different than package X" queries.  

  Asking for results that are too close to machine precision
  (double.eps) is a fool's errand; a reasonable critera is often the
  square root of that precision.  The Cholesky decompostion needs to be
  held to a higher standard than the overall convergence criterion, however.
  The \code{tolerance.inf} value controls a warning message; if it is
  too small incorrect warnings can appear, if too large some actual cases of
  an infinite coefficient will not be detected.

  The most difficult cases are data sets where the MLE coefficient is
  infinite; an example is a data set where at each death time,
  it was the subject with the largest
  covariate value who perished.  In that situation the coefficient
  increases at each iteration while the log-likelihood asymptotes to a
  maximum.  As iteration proceeds there is a race condition
  condition for three endpoint: exp(coef) overflows,
  the Hessian matrix become singular, or the change in loglik is small
  enough to satisfy the convergence criterion.  The first two are
  difficult to anticipate and lead to numeric diffculties, which is
  another argument for moderation in the choice of \code{eps}.

  See the vignette "Roundoff error and tied times" for a more
  detailed explanation of the \code{timefix} option.  In short, when
  time intervals are created via subtraction then two time intervals that are
  actually identical can appear to be different due to floating point
  round off error, which in turn can make \code{coxph} and
  \code{survfit} results dependent
  on things such as the order in which operations were done or the
  particular computer that they were run on.
  Such cases are unfortunatedly not rare in practice.
  The \code{timefix=TRUE} option adds
  logic similar to \code{all.equal} to ensure reliable results.
  In analysis of simulated data sets, however, where often by defintion there
  can be no duplicates, the option will often need to be set to 
  \code{FALSE} to avoid spurious merging of close numeric values.
}  
\seealso{\code{\link{coxph}}
  }
\keyword{survival}
