\name{survfit23}
\alias{survfit23}
\alias{survfit32}
\title{
Convert from old to new format for survfit objects
}
\description{
  In version 3.x of the package the components of a \code{survfit}
  object changed form.  These two routines translate back and forth
  between them.
}
\usage{
survfit23(x)
survfit32(x)
}
\arguments{
  \item{x}{a survfit object}
}
\details{
  A design decision when \code{survfit} objects were first created
  (approximately 1991) was to not include the first (x=0, y=1) data
  point of the curve in the \code{time} and \code{surv} portions of
  the results.  That turned out to be unwise.
  With extentions to survfit methods over the years,
  the first data point is not always (0,1), which neccessitated adding
  a sequence of additional components to hold that information,
  e.g., \code{start.time}, \code{p0}, and \code{std.p0}.
  Addition of the cumulative hazard and threatened to add even more, and
  in version 3.x of the package the object changed.

  In the newer version \code{start.time} is absorbed as the first
  element of the \code{time} vector, \code{p0} added as the first
  row of \code{pstate}, etc.  A flag \code{version=3} is also added
  to the object as a marker.   Since this transtion will not happen all
  at once for all component of the library, and for the use of other
  users, the \code{survfit23} routine can be used to convert from
  survival 2.x form to survival 3.x form, and vice-versa with
  \code{survfit32}.  If applied to an object that is already in the
  desired form they do nothing.
}
\value{ a survfit object}
\keyword{survival}
