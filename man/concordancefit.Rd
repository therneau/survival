\name{concordancefit}
\alias{concordancefit}
\title{Compute the concordance}
\description{
This is the working routine behind the \code{concordance} function.  It
is not meant to be called by users, but is available for other packages
to use.  Input arguments, for instance, are assumed to all be the
correct length and type, and missing values are not allowed: the calling
routine is responsible for these things.
}
\usage{
concordancefit(y, x, strata, weights, ymin = NULL, ymax = NULL,
 timewt = c("n", "S", "S/G", "n/G2", "I"), cluster, influence =0,
 ranks = FALSE, reverse = FALSE, timefix = TRUE, keepstrata=10, 
 std.err = TRUE)
}
\arguments{
  \item{y}{the response.  It can be numeric, factor, or a Surv object}
  \item{x}{the predictor, a numeric vector}
  \item{strata}{optional numeric vector that stratifies the data}
  \item{weights}{options vector of case weights}
  \item{ymin, ymax}{restrict the comparison to response values in this
    range}
  \item{timewt}{the time weighting to be used}
  \item{cluster, influence,ranks, reverse, timefix}{see the help for the
    \code{concordance} function}
  \item{keepstrata}{either TRUE, FALSE, or an integer value.
    Computations are always done within stratum, then added. If the
    total number of strata greater than \code{keepstrata}, or
    \code{keepstrata=FALSE}, those subtotals are not kept in the output.
  }
  \item{std.err}{compute the standard error; not doing so saves some
  compute time.}
}
\details{
  This function is provided for those who want a ``direct'' call to the
  concordance calculations, without using the formula interface.  A
  primary use has been other packages.   The routine does minimal
  checking of its input arguments, under the assumption that this has
  already been taken care of by the calling routine.
}
\value{a list containing the results}
\author{ Terry Therneau}
\seealso{ \code{\link{concordance}}}
\keyword{ survival }
