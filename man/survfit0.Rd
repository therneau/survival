\name{survfit0}
\alias{survfit0}
\title{
Convert the format of a survfit object.
}
\description{
  Add the point for a starting time (time 0) to a survfit object's
  elements.  This is useful for plotting.
}
\usage{
survfit0(x, start.time=0)
}
\arguments{
  \item{x}{a survfit object}
  \item{start.time}{the desired starting time}
}
\value{a reformulated version of the object with an initial data point
  at \code{start.time} added. 
  The \code{time}, \code{surv}, \code{pstate}, \code{cumhaz},
  \code{std.err}, and \code{std.cumhaz} components will all be aliged,
  so as to make plots and summaries easier to produce.
}

\details{
  Survival curves are traditionally plotted forward from time 0, but
  since the true starting time is not known as a part of the data,
  the \code{survfit} routine does not include a time 0 value in
  the resulting object.
  Someone might look at cumulative mortgage defaults versus calendar
  year, for instance, with the `time' value a Date object.
  The plotted curve probably should not start at 0 = 1970/01/01.
  (One could certainly argue that this behavior is overly pedantic, but
  the decision was made over 20 years ago and it is far too late to
  change it now.)
  
  This is a utility function used by \code{plot.survfit} and
  \code{summary.survfit}.  If the original \code{survfit} call
  included a \code{start.time} argument then that will be used as
  the starting point, otherwise the default of 0.
  It may be useful to others, so was exported.
}
\keyword{survival}
