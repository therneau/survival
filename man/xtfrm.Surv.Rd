\name{xtfrm.Surv}
\alias{xtfrm.Surv}
\alias{sort.Surv}
\alias{order.Surv}
\title{Sorting order for Surv objects}
\description{
  Sort survival objects into a partial order, which is the same one
  used internally for many of the calculations.  
}
\usage{
\method{xtfrm}{Surv}(x)
}
\arguments{
  \item{x}{a \code{Surv} object}
}
\details{
  This creates a partial ordering of survival objects.
  The result is sorted in time order, for tied pairs of times right censored
  events come after observed events (censor after death), and left
  censored events are sorted before observed events.
  For counting process data \code{(tstart, tstop, status)} the ordering
  is by stop time, status, and start time, again with censoring last.
  Interval censored data is sorted using the midpoint of each interval.
  
  The \code{xtfrm} routine is used internally by \code{order} and
  \code{sort}, so these results carry over to those routines.
}
\value{a vector of integers which will have the same sort order as
  \code{x}.
}
\author{Terry Therneau}
\seealso{\code{\link{sort}}, \code{\link{order}}}
\examples{
test <- c(Surv(c(10, 9,9, 8,8,8,7,5,5,4), rep(1:0, 5)), Surv(6.2, NA))
test
sort(test)
}
\keyword{survival}

