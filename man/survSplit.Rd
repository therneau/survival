\name{survSplit}
\alias{survSplit}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Split a survival data set at specified times }
\description{
  Given a survival data set and a set of specified cut times, split
  each record into multiple subrecords at each cut time.  The new data
  set will be in `counting process' format, with a start time, stop
  time, and event status for each record. 
}
\usage{
survSplit(data, cut, end, event, start, id = NULL, zero = 0,
  episode=NULL)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{data}{data frame}
  \item{cut}{vector of timepoints to cut at}
  \item{end}{character string with name of event time variable }
  \item{event}{character string with name of censoring indicator }
  \item{start}{character string with name of start time variable (will
    be created if it does not exist) }
  \item{id}{character string with name of new id variable to create (optional)}
  \item{zero}{If \code{start} doesn't already exist, this is the time
    that the original records start. May be a vector or single value.}
  \item{episode}{character string with name of new episode variable (optional)}
}

\value{
 New, longer, data frame.
}

\details{
  The function also works when the original data are in
  counting-process format, but the \code{id} and \code{episode} options
  are of little use in this context.
 
  }

\seealso{\code{\link{Surv}}, \code{\link{cut}}, \code{\link{reshape}} }

\examples{
aml3<-survSplit(aml,cut=c(5,10,50),end="time",start="start",
      event="status",episode="i")

summary(aml)
summary(aml3)

coxph(Surv(time,status)~x,data=aml)
## the same
coxph(Surv(start,time,status)~x,data=aml3)

aml4<-survSplit(aml3,cut=20,end="time",start="start", event="status")
coxph(Surv(start,time,status)~x,data=aml4)

}
\keyword{survival }% at least one, from doc/KEYWORDS
\keyword{utilities}% __ONLY ONE__ keyword per line
