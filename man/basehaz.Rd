\name{basehaz}
\alias{basehaz}
\title{Alias for the survfit function}
\description{
  Compute the predicted survival curve for a Cox model.
}
\usage{
basehaz(\dots)
}
\arguments{
  \item{\dots}{arguments to \code{survfit}}
}
\details{
  This function is simply an alias for \code{survfit}, which is the
  actual function that does all the computations.  See the manual page
  for that function for the remainder of the arguments.
  Since survival curves after a Cox model are based on the baseline
  hazard function, users often look for a function of that name.

  Earlier versions of this function would return only selected portions
  of the survfit object, and also renamed the cumulative hazard portion as
  'hazard'.
}
\value{
  the result of a call to \code{survfit}
}

\seealso{\code{\link{survfit}}}
\keyword{survival }

