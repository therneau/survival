\name{survfit0}
\alias{survfit0}
\title{
Convert the format of a survfit object.
}
\description{
  Add the point for a starting time ("time 0") to a survfit object's
  elements.  This is useful for plotting.
}
\usage{
survfit0(x, ...)
}
\arguments{
  \item{x}{a survfit object}
  \item{...}{any other arguments are ignored}
  }
\value{a reformulated version of the object with an initial data point added.
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
  The plotted curve probably should not start at 0 = 1970-01-01.
  Due to this uncertainty, it was decided not to include a "time 0" as
  part of a survfit object. 
  Whether that (1989) decision was wise or foolish,
  it is now far too late to change it. (We tried it once as a
  trial, resulting in over 20 errors in the survival test suite.  We
  extrapolated that it might break 1/3 of the other CRAN packages
  that depend on survival, if made a default.)
  Many curves do include a value \code{t0} for "time 0",
  which is where the survfit
  routine has surmised that the curve would start.
    
  One problem with this choice is that some functions must choose a
  starting point, plots and computation of the restricted mean survival
  time are two primary examples.
  This utility function is used by \code{plot.survfit} and
  \code{summary.survfit} to fill in that gap.

  The value used for this first time point is the first one below
  \enumerate{
    \item a \code{t0} value found in the  in the object. 
    \item for single state survival \itemize{
      \item min(0, time) for Surv(time, status) data
      \item min(time1) for Surv(time1, time2, status) data
    }
    \item for multi state survival \itemize{
      \item min(0, time) for Surv(time, event) data, e.g., competing
      risks
      \item min(time1) for Surv(time1, time2, event) data, if everyone
      starts in the same state
    }
  }
  (Remember that negative times are allowed in Surv objects.)
  
  This function will add a new time point at the front of each curve,
  but only if said time point is less than existing points in the
  curve.  If there were a death on day 0, for instance, it will not add a
  (time=0, survival=1) point.
  (The question of whether the plotted curve in this case should or
  should not start with a vertical segment can be debated ad nauseum.
  It has no effect on the area under the curve (RMST), and the summary
  for time 0 should report the smaller value.) 

  The resulting object is \emph{not} currently
  guarranteed to work with functions that further manipulate a
  \code{survfit} object such as subscripting, aggregation, pseudovalues,
  etc. (remember the 20 errors).  Rather it is intended as a penultimate
  step, most often when creating a plot or summary of the curve(s).
}
\keyword{survival}
