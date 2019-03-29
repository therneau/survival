\name{multicheck}
\alias{multicheck}
\title{Checks of a survival data set}
\description{
Perform a set of consistency checks on survival data
}
\usage{
multicheck(y, id, istate = NULL, dummy = "( )")
}
\arguments{
  \item{y}{a \code{Surv} object}

  \item{id}{vecor of subject identifiers}
  \item{istate}{initial state for each observation, can be missing}
  \item{dummy}{the label that will be used for the "unnamed state"
    in the result.}
}
\details{
  This is called internally by several routines.
  It detects several conditions.
  \describe{
    \item{overlap}{two observations for the same subject that overlap in
      time, e.g. intervals of (0, 100) and (90, 120).  A subject cannot
      be two places at once.
      If \code{y} is simple (time, status) survival, observation
      intervals implicitly start at 0, so any duplicate identifiers will
      generate and overlap.}
    \item{gaps}{gaps in a subject's timeline, with no \code{istate}
      variable to resolve issues.  (Where were they?)}
    \item{inconsistent}{a subject who transitioned to one state on a
      particular day, based on the \code{y} variable, but appears in a
      different state, on the same day, based in \code{istate}.}
      \item{jumps}{a subject with holes in their follow-up.}
    }
    The total number of occurences is present in the \code{flags}
  vector. Optional components give the location and identifiers of the
  flagged observations.
}
\value{
  a list with components
  \item{transitions}{a matrix giving the count of transitions from one
    state to another}
  \item{istate}{a revised initial state vector}
  \item{states}{vector of state names}
  \item{flags}{a vector giving the counts of errors}
  \item{overlap}{a list with the row number and id of overlaps}
  \item{gaps}{a list with the row number and id of gaps}
  \item{inconsistent}{a list with the row number and id of inconsistent
    rows}
  \item{jumps}{a list with the row number and id of jumps}
}
\author{Terry Therneau}

\keyword{ survival }

