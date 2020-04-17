\name{Surv2data}
\alias{Surv2data}
\title{Convert data from timecourse to (time1,time2) style
}
\description{
  The multi-state survival functions \code{coxph} and \code{survfit}
  allow for two forms of input data.  This routine converts between them.
  The function is normally called behind the scenes when \code{Surv2} is
  as the response. 
}
\usage{
Surv2data(formula, data, subset, id)
}
\arguments{
  \item{formula}{a model formula}
  \item{data}{a data frame}
  \item{subset}{optional, selects rows of the data to be retained}
  \item{id}{a variable that identified multiple rows for the same
    subject, normally found in the referenced data set}
}
\value{
  a list with elements
  \item{mf}{an updated model frame (fewer rows, unchanged columns)}
  \item{S2.y}{the constructed response variable}
  \item{S2.state}{the current state for each of the rows}
}
  
\details{
For timeline style data, each row is uniquely identified by an
(identifier, time) pair.  The time could be a date, time from entry to a
study, age, etc, (there may often be more than one time variable).
The identifier and time cannot be missing.
The remaining covariates represent values that were observed at that
time point. Often, a given covariate is observed at only a subset of
times and is missing at others.  At the time of death, in particular,
often only the identifier, time, and status indicator are known.

In the resulting data set missing covariates are replaced by their
last known value, and the response y will be a Surv(time1, time2,
endpoint) object.
}
\keyword{survival}
