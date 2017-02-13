\name{print.summary.survexp}
\alias{print.summary.survexp}
\title{Print Survexp Summary}
\description{
Prints the results of \code{summary.survexp}
}
\usage{
\method{print}{summary.survexp}(x, digits = max(options()$digits - 4, 3), ...)
}
\arguments{
  \item{x}{
an object of class \code{summary.survexp}.
}
  \item{digits}{
the number of digits to use in printing the result.
}
  \item{\dots}{for future methods}
}
\value{
 \code{x}, with the invisible flag set to prevent further printing.
}
\author{Terry Therneau}
\seealso{\code{link{summary.survexp}}, \code{\link{survexp}}}
\keyword{ survival }

