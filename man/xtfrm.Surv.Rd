\name{xtfrm.Surv}
\alias{xtfrm.Surv}
\alias{sort.Surv}
\title{Sorting order for Surv objects}
\description{
  Sort survival objects into a partial order, which is the same one
  used internally for many calculations.  
}
\usage{
\method{xtfrm}{Surv}(x)
\method{sort}{Surv}(x, \ldots)
}
\arguments{
  \item{x}{a \code{Surv} object}
  \item{\ldots}{other arguments to the \code{sort} function}
}
\details{
  This creates a partial ordering of survival objects.
  The results is sorted in time order, for tied times censored
  events come after observed events (censor after death).
  For counting process data \code{(tstart, tstop, status)} the ordering
  is by stop time, start time, status, again with censoring last.
  The \code{xtfrm} routine is used internally by \code{order}.
}
\value{a vector of integers (\code{xtfrm}) or the ordered survival times
  (\code{sort}).
}
\author{Terry Therneau}
\seealso{\code{\link{sort}}}
\examples{
test <- c(Surv(c(10, 9,9, 8,8,8,7,5,5,4), rep(1:0, 5)), Surv(6.2, NA))
test
sort(test)
}
\keyword{survival}

