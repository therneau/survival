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
survSplit(data, cut, end, event, start="tstart", id, zero = 0,
  episode)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{data}{a data frame}
  \item{cut}{the vector of timepoints to cut at}
  \item{end}{character string with the name of event time variable }
  \item{event}{character string with the name of censoring indicator }
  \item{start}{character string with the name of the start time variable (will
    be created if it does not exist) }
  \item{id}{character string with the name of new id variable to
    create (optional).  This can be useful if the data set does not
    already contain an identifier.}
  \item{zero}{If \code{start} doesn't already exist, this is the time
    that the original records start. May be a vector or single value.}
  \item{episode}{character string with the name of new episode
    variable (optional)}
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
  }

\seealso{\code{\link{Surv}}, \code{\link{cut}}, \code{\link{reshape}} }

\examples{
fit1 <- coxph(Surv(time, status) ~ karno + age + trt, veteran)
# a cox.zph plot of the data suggests that the effect of Karnofsky score
#  has faded by 120 days-- plot(cox.zph(fit)[1]) -- and begins to
#  diminish by 60 days.
# Fit a model with separate coefficeints for the three intervals
#
vet2 <- survSplit(veteran, cut=c(60, 120), end="time", event="status",
                  episode ="timegroup")
fit2 <- coxph(Surv(tstart, time, status) ~ karno* strata(timegroup) +
                age + trt, data= vet2)
c(overall= coef(fit1)[1],
  t0_60  = coef(fit2)[1],
  t60_120= sum(coef(fit2)[c(1,4)]),
  t120   = sum(coef(fit2)[c(1,5)]))

}
\keyword{survival }
\keyword{utilities}
