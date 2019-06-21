\name{survcheck}
\alias{survcheck}
\title{Checks of a survival data set}
\description{
Perform a set of consistency checks on survival data
}
\usage{
survcheck(formula, data, id, istate, istate0="(s0)", ...)
}
\arguments{
  \item{formula}{a model formula with a \code{Surv} object as the
    response}
  \item{data}{data frame in which to find the \code{id},
    \code{istate} and formula variables}
  \item{id}{an identifier that labels unique subjects}
  \item{istate}{an optional vector giving the initial state for each
    observation, i.e., the state at the start of each interval}
  \item{istate0}{default label for the initial state when \code{istate}
    is missing}
  \item{\ldots}{other arguments, which are ignored (but won't give an
    error if someone added \code{weights} for instance)}
}
\details{
  This routine will examine a multi-state data set for consistency of
  the data.  The basic rules are that if a subject is at risk they have
  to be somewhere, can not be two states at once, and should make
  sensible transitions from state to state.  It reports the number of
  instances of the following conditions:
  \describe{
    \item{overlap}{two observations for the same subject that overlap in
      time, e.g. intervals of (0, 100) and (90, 120).  
      If \code{y} is simple (time, status) survival observation
      intervals implicitly start at 0, so in that case any duplicate
      identifiers will generate an overlap.}
    \item{jump}{a hole in a subject's timeline, where they are in one
      state at the end of the prior interval, but a new state in the
      at the start subsequent interval.}
   \item{gap}{one or more gaps in a subject's timeline; they are presumably
      in the same state at their return as when they left.}
    \item{teleport}{two adjacent intervals for a subject, with the
      first interval ending in one state and the subsequent interval
      starting in another.  They have instantaneously changed states
      with experiencing a transition.}
     }
  The total number of occurences of each is present in the \code{flags}
  vector. Optional components give the location and identifiers of the
  flagged observations.
}
\value{
  a list with components
  \item{states}{the vector of possible states}
  \item{transitions}{a matrix giving the count of transitions from one
    state to another}
  \item{statecount}{table of the number of visits per state, e.g., 18
    subjects had 2 visits to the "infection" state}
  \item{flags}{a vector giving the counts of each check}
  \item{istate}{a revised initial state vector}
  \item{overlap}{a list with the row number and id of overlaps (not
    present if there are no overlaps)}
  \item{gaps}{a list with the row number and id of gaps (not present if
    there are no gaps)}
  \item{teleport}{a list with the row number and id of inconsistent
    rows (not present if there are none)}
  \item{jumps}{a list with the row number and id of jumps (not present
    if there are no jumps}
}
\author{Terry Therneau and Beth Atkinson}

\keyword{ survival }

