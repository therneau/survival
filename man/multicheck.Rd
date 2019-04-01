\name{multicheck}
\alias{multicheck}
\title{Checks of a survival data set}
\description{
Perform a set of consistency checks on survival data
}
\usage{
multicheck(formula, data, id, istate, ...)
}
\arguments{
  \item{formula}{a model formula with a \code{Surv} object as the
    response}
  \item{data}{data frame in which to find the \code{id},
    \code{istate} and formula variables}
  \item{id}{an identifier that labels unique subjects}
  \item{istate}{an optional vector giving the initial state for each
    observation, i.e., the state at the start of each interval}
  \item{\ldots}{other arguments, which are ignored (but don't give an
    error if someone added \code{weights} for instance)}
}
\details{
  This routine will examine a multi-state data set for consistency of
  the data.  The basic rules are that if a subject is at risk they have
  to be somewhere, but can not be two places at once.
  
  It detects several conditions.
  \describe{
    \item{overlap}{two observations for the same subject that overlap in
      time, e.g. intervals of (0, 100) and (90, 120).  A subject cannot
      be two places at once.
      If \code{y} is simple (time, status) survival, observation
      intervals implicitly start at 0, so any duplicate identifiers will
      generate and overlap.}
    \item{gap}{one or more gaps in a subject's timeline, with no \code{istate}
      variable to resolve issues.  (Where were they?)}
    \item{teleport}{a subject who transitioned to one state on a
      particular day, based on the \code{y} variable, but appears in a
      different state, on the same day, based in \code{istate}.}
    \item{jump}{a subject with holes in their follow-up.}
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

