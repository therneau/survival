\name{aeqSurv}
\alias{aeqSurv}
\title{Adjudicate near ties in a Surv object}
\description{The check for tied survival times can fail due
  to floating point imprecision, which can make actual ties appear to
  be distinct values.
  Routines that depend on correct identification of ties pairs will then
  give incorrect results, e.g., a Cox model.
  This function rectifies these.
}
\usage{
aeqSurv(x, tolerance = sqrt(.Machine$double.eps))
}
\arguments{
  \item{x}{a Surv object}
  \item{tolerance}{the tolerance used to detect values that will
    be considered equal}
}
\details{
  This routine is called by both \code{survfit} and \code{coxph} to
  deal with the issue of ties that get incorrectly broken due to
  floating point imprecision.  See the short vignette on this issue
  for a simple example.  Use the \code{timefix} argument of
  \code{survfit} or \code{\coxph.control} to control the option
  if desired.

  The rule for `equality' is identical to that used by the
  \code{all.equal} routine.  Pairs of values that are within round off
  error of each other are replaced by the smaller value.
  An error message is generated if this process causes a 0 length
  time interval to be created.  
}
\value{a Surv object identical to the original, but with ties restored.}
\author{Terry Therneau}
\seealso{\code{\link{survfit}}, \code{\link{coxph.control}}}

\keyword{ survival }

