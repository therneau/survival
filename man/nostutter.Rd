\name{nostutter}
\alias{nostutter}
\title{
Remove repeated events
}
\description{
For a subject with repeated follow-up visits, do not count the second,
third, etc. visit in a row with the same state as a repeated event.
}
\usage{
nostutter(id, x, censor = 0)
}
\arguments{
  \item{id}{subject identifier}
  \item{x}{covariate of interest}
  \item{censor}{censoring code to be used as the "no event at this time"
    value}
}
\details{
  When an covariate is recorded multiple times during follow-up, and that
  variable is being used as an endpoint in a multi-state hazard model,
  we will sometimes want to consider each instance as a new event and
  sometimes not. An example of the first kind is the \code{cgd} data
  set; each row in which \code{status =1} is a new infection.
  An example of the second is the \code{nafld} data when we treat the
  number of co-morbidities as a state; multiple visits with the same
  total do not count as entry into a new state.

  This routine will often be called on a current state variable before
  using that variable as an event variable in \code{tmerge}.
}
\value{an updated vector with repeated values marked as \code{censor}}
\seealso{\code{\link{lvcf}}, \code{\link{tmerge}}
}
\keyword{ survival }
