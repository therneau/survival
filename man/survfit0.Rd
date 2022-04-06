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
  \item{start.time}{the desired starting time; see details below.}
  }
\value{a reformulated version of the object with an initial data point
  at \code{start.time} added. 
  The \code{time}, \code{surv}, \code{pstate}, \code{cumhaz},
  \code{std.err}, \code{std.cumhaz} and other components will all be aligned,
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
  Due to this uncertainty, it was decided not to include a "time 0" as
  part of a survfit object.  If the original \code{survfit} call
  included a \code{start.time} argument, that value is of course
  retained.

  Whether that (1989) decision was wise or
  foolish, it is now far too late to change it. (We tried it once as a
  trial, resulting in over 20 errors in the survival test suite.  We
  extrapolate that it might break 1/2 - 2/3 of the other CRAN packages
  that depend on survival, if made a default.)
  If the original \code{survfit} call
  included a \code{start.time} argument, that value is of course
  retained.
    
  One problem with this choice is that some functions must choose a
  starting point, plots for example.
  This utility function is used by \code{plot.survfit} and
  \code{summary.survfit} to do so, adding a new time point at the front of each
  curve in a consistent way: the optional argument to the \code{survfit0}
  function as the first choice (if supplied),
  then the user's \code{start.time} if present,
  otherwise \code{min(0, x$time)}.  The resulting object is \emph{not}
  guarranteed to work with functions that further manipulate a
  \code{survfit} object such as subscripting, aggregation, pseudovalues,
  etc. (remember the 20 errors).  Rather it is intended as a penultimate
  step, most often when creating a plot.
}
\keyword{survival}
