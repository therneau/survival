\name{survcondense}
\alias{survcondense}
\title{Shorten a (time1, time2) survival dataset
}
\description{
  Counting process data sets can sometimes grow to an uweildy size, this can
  be used to reduce the number of rows.
}
\usage{
survcondense(formula, data, subset, weights, na.action= na.pass, id, 
             start = "tstart", end = "tstop", event = "event")
}
\arguments{
  \item{formula}{
    a formula object, with the response on the left of a \code{~} operator, and 
    the terms on the right.  The response must be a survival object as 
    returned by the \code{Surv} function. 
  }
  \item{data}{
    a data.frame in which to interpret the variables named in 
    the \code{formula} and the \code{id} argument
    argument. 
  }
  \item{subset}{optional subset expression to apply to the data set}
  \item{weights}{optional variable name for case weights}
  \item{na.action}{optional removal of missing values}
  \item{id}{variable name that identifies subjects}
  \item{start}{optional character string, giving the name of the start
    time variable in the result}
  \item{end}{optional character string, giving the name of the stop
    time variable in the result}
  \item{event}{optional character string, giving the name of the event
    variable in the result}
}
\details{
  Through the use of the \code{survSplit} and \code{tmerge} functions, a
  counting process data set will gain more and more rows of data.
  Occassionally it is useful to collapse this surplus back down, e.g.,
  when interest is to be focused on only a few covariates, or for
  debugging.  The right hand side of \code{formula} will often have only
  a few variables in this use case.

  If a row of data is censored, and represents the same covariates and
  identifier as the row below it, then the two rows can be merged
  together using a single (time1, time2) interval.  The compression can
  sometimes be large.

  The \code{start}, \code{stop} and \code{end} options are only used when the
  left hand side of the formula has expressions that are not a simple
  name, e.g. \code{Surv(time1, time2, death | progression)} would be a
  case where \code{event} is used to set the outcome variable's name.
}
\value{a data frame}
\author{Terry Therneau}
\seealso{
\code{\link{survSplit}},\code{\link{tmerge}}
}
\examples{
dim(aml)
test1 <- survSplit(Surv(time, status) ~ ., data=aml, 
                   cut=c(10, 20, 30), id="newid")
dim(test1)

# remove the added rows
test2 <- survcondense(Surv(tstart, time, status) ~ x, test1, id=newid)
dim(test2)
}
\keyword{ survival }
