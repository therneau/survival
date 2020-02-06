\name{survcheck}
\alias{survcheck}
\title{Checks of a survival data set}
\description{
Perform a set of consistency checks on survival data
}
\usage{
survcheck(formula, data, id, istate, istate0="(s0)", timefix=TRUE,...)
}
\arguments{
  \item{formula}{a model formula with a \code{Surv} object as the
    response}
  \item{data}{data frame in which to find the \code{id},
    \code{istate} and formula variables}
  \item{id}{an identifier that labels unique subjects}
  \item{istate}{an optional vector giving the current state at the start
    of each interval}
  \item{istate0}{default label for the initial state of each subject (at
    their first interval) when \code{istate} is missing}
  \item{timefix}{process times through the \code{aeqSurv} function to
  eliminate potential roundoff issues.}
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
  \item{istate}{a copy of the istate vector, if it was supplied;
    otherwise a constructed istate that satisfies all the checks}
  \item{overlap}{a list with the row number and id of overlaps (not
    present if there are no overlaps)}
  \item{gaps}{a list with the row number and id of gaps (not present if
    there are no gaps)}
  \item{teleport}{a list with the row number and id of inconsistent
    rows (not present if there are none)}
  \item{jumps}{a list with the row number and id of jumps (not present
    if there are no jumps}
}
\note{
For data sets with time-dependent covariates, a given subject will often
have intermediate rows with a status of 0 (no event at this time), for
instance a subject who started in state 1 at time 0, transitioned to state
2 at time 10, had a covariate \code{x} change from 135 to 156 at time
20, and a final transition to state 3 at time 30.
The response would be \code{Surv(c(0, 10, 20), c(10, 20, 30), c(1,0,2))}:
the status variable records \emph{changes} in state, and there was no
change at time 20.
The \code{istate} variable would be (1, 2, 2); it retains its prior
value status = censored.

Thus, when there are intermediate observations \code{istate} is not
simply a lagged version of the status, and may be more challenging to
create.
A simple approach is to call the \code{survcheck} routine without an
\code{istate} argument, and use the constructed \code{istate} variable
that it returns.  

\keyword{ survival }

