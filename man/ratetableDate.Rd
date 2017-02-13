\name{ratetableDate}
\alias{ratetableDate}
\title{Convert date objects to ratetable form}
\description{
  This method converts dates from various forms into
  the internal form used in \code{ratetable} objects.
}
\usage{
ratetableDate(x)
}
\arguments{
  \item{x}{a date.  The function currently has methods for Date, date,
    POSIXt, timeDate, and chron objects.
}
}
\details{
  This function is useful for those who create new ratetables, but is
  normally invisible to users.
  It is used internally by the \code{survexp} and \code{pyears}
  functions to map the various date formats; if a new method is added
  then those routines will automatically be adapted to the new date type.
}
\value{a numeric vector, the number of days since 1/1/1960.}
\author{Terry Therneau}
\seealso{\code{\link{pyears}}, \code{\link{survexp}}}
\keyword{survival}

