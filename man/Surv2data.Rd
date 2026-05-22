\name{Surv2data}
\alias{Surv2data}
\title{Convert data from timecourse to (time1,time2) style
}
\description{
  The multi-state survival functions \code{coxph} and \code{survfit}
  allow for two forms of input data.  This routine converts between them.
  The function is normally called behind the scenes when \code{Surv2} is
  used as the response. 
}
\usage{
Surv2data(formula, data, subset, id, lcvf=TRUE)
}
\arguments{
  \item{formula}{a model formula}
  \item{data}{a data frame}
  \item{subset}{optional, selects rows of the data to be retained}
  \item{id}{a variable that identifies multiple rows for the same
    subject, necessary if there are multiple rows per subject}
}
\value{ a model frame}  
\details{
For timeline style data, each row is uniquely identified by an
(identifier, time) pair.  The time could be a date, time from entry to a
study, age, etc.
The identifier and time variables cannot be absent or have missing values.
The remaining covariates represent values that were observed at that
time point. Often, a given covariate is observed at only a subset of
times and is missing at others.  At the time of death, in particular,
often only the identifier, time, and status indicator are known.
Covariates such as sex or genetic markers, which are constant over time,
might be present in all rows.

The survfit and coxph functions allow data in either
counting process form (the historic format) or timeline form; the user
is responsible for using \code{Surv2(time, status)} to identify that this
is timeline data, as opposed to counting process
(\code{Surv(time1, time2, status)}) or simple (\code{time, status}) survival
data.
Internally, they immediately convert the data to standard or
counting process form as appropriate. This new data set does not have
the same number of rows.
The results of the fit correspond to that transformed data set, e.g.,
the \code{linear.predictors} and \code{model} components of the fit,
or predicted values and residuals from the \code{predict} and
\code{residuals} functions.

This is a utility routine, it returns a model frame in the updated form,
containing an extra variable \code{_orow_}, the row number in the orignal
data set.
\keyword{survival}
