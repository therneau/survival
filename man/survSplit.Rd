\name{survSplit}
\alias{survSplit}
\title{Split a survival data set at specified times }
\description{
  Given a survival data set and a set of specified cut times, split
  each record into multiple subrecords at each cut time.  The new data
  set will be in `counting process' format, with a start time, stop
  time, and event status for each record. 
}
\usage{
survSplit(formula, data, subset, na.action=na.pass,
            cut, start="tstart", id, zero=0, episode,
                              end="tstop", event="event")
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{formula}{a model formula}
  \item{data}{a data frame}
  \item{subset, na.action}{rows of the data to be retained}
  \item{cut}{the vector of timepoints to cut at}
  \item{start}{character string with the name of a start time variable (will
    be created if needed) }
  \item{id}{character string with the name of new id variable to
    create (optional).  This can be useful if the data set does not
    already contain an identifier.}
  \item{zero}{If \code{start} doesn't already exist, this is the time
    that the original records start.}
  \item{episode}{character string with the name of new episode
    variable (optional)}
  \item{end}{character string with the name of event time variable }
  \item{event}{character string with the name of censoring indicator }
}

\value{
 New, longer, data frame.
}

\details{
  Each interval in the original data is cut at the given points; if
  an original row were (15, 60] with a cut vector of (10,30, 40) the
  resulting data set would have intervals of (15,30], (30,40] and
  (40, 60].

  Each row in the final data set will lie completely within one of the
  cut intervals. Which interval for each row of the output is shown by the
  \code{episode} variable, where 1= less than the first cutpoint, 2=
  between the first and the second, etc.
  For the example above the values would be 2, 3, and 4.

  The routine will normally be called with a formula as the first
  argument. 
  The right hand side of the formula can be used to delimit variables
  that should be retained; normally one will use \code{ ~ .} as a
  shorthand to retain them all.  When called in this way the routine
  will try to retain variable names, e.g. \code{Surv(adam, joe, fred)~.}
  will result in a data set with those same variable names for
  \code{tstart}, \code{end}, and \code{event} options rather than
  the defaults.  Any user specified values for these options will be
  used if they are present, of course.
  However, the routine is not sophisticated; it only does this
  substitution for simple names.  A call of \code{Surv(time, stat==2)}
  for instance will result in use of the 'event' argument for that
  variable name.
  The \code{end} and \code{event} arguments are required if there is
  no formula.

  Rows of data with a missing time or status are copied across
  unchanged, unless the na.action argument is changed from its default
  value of \code{na.pass}.  But in the latter case any row
  that is missing for any variable will be removed, which is rarely
  what is desired.
  }

\seealso{\code{\link{Surv}}, \code{\link{cut}}, \code{\link{reshape}} }

\examples{
fit1 <- coxph(Surv(time, status) ~ karno + age + trt, veteran)
plot(cox.zph(fit1)[1])
# a cox.zph plot of the data suggests that the effect of Karnofsky score
#  begins to diminish by 60 days and has faded away by 120 days.
# Fit a model with separate coefficients for the three intervals.
#
vet2 <- survSplit(Surv(time, status) ~., veteran,
                   cut=c(60, 120), episode ="timegroup")
fit2 <- coxph(Surv(tstart, time, status) ~ karno* strata(timegroup) +
                age + trt, data= vet2)
c(overall= coef(fit1)[1],
  t0_60  = coef(fit2)[1],
  t60_120= sum(coef(fit2)[c(1,4)]),
  t120   = sum(coef(fit2)[c(1,5)]))
}

\keyword{survival }
\keyword{utilities}
