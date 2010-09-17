\name{tcut}
\alias{tcut}
\alias{[.tcut}
\title{Factors for person-year calculations}
\description{
Attaches categories for person-year calculations to a variable without
losing the underlying continuous representation
}
\usage{
tcut(x, breaks, labels, scale=1)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{numeric/date variable }
  \item{breaks}{breaks between categories, which are right-continuous }
  \item{labels}{labels for categories }
  \item{scale}{Multiply \code{x} and \code{breaks} by this.}
}

\value{
  An object of class \code{tcut}
}

\seealso{ \code{\link{cut}}, \code{\link{pyears}} }

\examples{
mdy.date <- function(m,d,y)
   as.Date(paste(ifelse(y<100, y+1900, y), m, d, sep='/'))
temp1 <- mdy.date(6,6,36)
temp2 <- mdy.date(6,6,55)# Now compare the results from person-years
#
temp.age <- tcut(temp2-temp1, floor(c(-1, (18:31 * 365.24))),
	labels=c('0-18', paste(18:30, 19:31, sep='-')))
temp.yr  <- tcut(temp2, mdy.date(1,1,1954:1965), labels=1954:1964)
temp.time <- 3700   #total days of fu
py1 <- pyears(temp.time ~ temp.age + temp.yr, scale=1) #output in days
py1
}
\keyword{survival}

