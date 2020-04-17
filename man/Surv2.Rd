\name{Surv2}
\alias{Surv2}
\title{Create a survival object}
\description{
  Create a survival object from a timeline style data set. This will
  almost always be the response variable in a formula.
}
\usage{
Surv2(time, event)
}
\arguments{
  \item{time}{a timeline variable, such as age, time from enrollment,
  date, etc.}
\item{event}{the outcome at that time.  This can be a 0/1 variable,
  TRUE/FALSE, or a factor.
  If the latter, the first level of the factor corresponds to
  `no event was observed at this time'.}
}

\value{
  An object of class \code{Surv2}.  There are methods for \code{print},
  \code{is.na} and subscripting.
}

\details{
  When used in a \code{coxph} or \code{survfit} model,
  Surv2 acts as a trigger to internally convert a timeline style data
  set into counting process style data, which is then acted on by the
  routine.
}
\seealso{
\code{\link{coxph}},  
\code{\link{survfit}}  
}
\keyword{survival}
