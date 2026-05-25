\name{timeline2counting}
\alias{timeline2counting}
\title{
Create counting process data from timeline data
}
\description{
The primary survival functions accept any of simple, counting process,
or timeline style arrangement of the data.  This converts the third to
the second
}
\usage{
timeline2counting(formula, data, subset, id, repeated = FALSE, lvcf = TRUE, 
                  yname = c("tstart", "tstop", "status"))
}
\arguments{
  \item{formula}{a formula with \code{Surv(time-variable,
      status-variable)} as the response. }
  \item{data}{a data set currently in timeline form}
  \item{subset}{optional expression to select a subset of rows}
  \item{id}{the identifier variable, which rows go with each subject}
  \item{repeated}{if TRUE, then two adjacent rows of data with the same
    response are treated as new events, repeated infections is one
    example.}
  \item{lvcf}{the last-value-carried-forward operation is applied to
    covariates on the right hand side}
  \item{yname}{names to be used for the new response variable set,
  see details below}
}
\details{
  The new data set will have fewer rows than the old, the row.names
  can be used to map back to the original.
  A new variable \code{istate} will be added, or \code{_istate_} if
  the former conflicts with a current variable.

  If the \code{yname} vector is supplied then those names will be used
  in the output.  If not, and the response has a simple structure,
  e.g. \code{Surv(age, state)} then the new data set will have variables
  of \code{age1, age2, state}, and the response when modeling with
  the new set would be \code{Surv(age1, age2, state)}.
  If the original response were more
  complex, \code{Surv(day, ifelse(day==0, 0, state2))} for instance,
  then the third element of \code{yname} would be used for the status
  variable. 
  If that that would generate a name conflict, i.e., the data set
  already had a variable named "status", then "_status_" would be used as
  the name.

  See the vignette on timeline data for a more complete discussion.
}
\value{a data frame}
\author{Terry Therneau}
\keyword{ survival }
