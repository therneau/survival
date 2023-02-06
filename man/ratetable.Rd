\name{ratetable}
\alias{ratetable}
\title{Allow ratetable() terms in a model}
\description{This function supports ratetable() terms in a model
  statement, within survexp and pyears.
}
\usage{
ratetable(...)
}
\arguments{
  \item{\dots}{the named dimensions of a rate table}
}
\details{
  This way of mapping a rate table's variable names to a user data frame
  has been superseded, instead use the \code{rmap} argument of the
  survexp, pyears, or survdiff routines.  The function remains only to
  allow older code to be run.
}
\author{Terry Therneau}
\keyword{survival}
