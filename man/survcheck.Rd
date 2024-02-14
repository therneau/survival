\name{survcheck}
\alias{survcheck}
\title{Checks of a survival data set}
\description{
Perform a set of consistency checks on survival data
}
\usage{
survcheck(formula, data, subset, na.action, id, istate, istate0="(s0)", 
timefix=TRUE,...)
}
\arguments{
  \item{formula}{a model formula with a \code{Surv} object as the
    response}
  \item{data}{data frame in which to find the \code{id},
    \code{istate} and formula variables}
  \item{subset}{expression indicating which subset of the rows of data
    should be used in the fit.    All observations are included by default. 
  }
  \item{na.action}{
    a missing-data filter function.  This is applied to the model.frame
    after any 
    subset argument has been used.  Default is \code{options()\$na.action}. 
  }
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
  to be somewhere, can not be two places at once, and should make
  sensible transitions from state to state. It reports the number of
  instances of the following conditions:
  \describe{
    \item{overlap}{two observations for the same subject that overlap in
      time, e.g. intervals of (0, 100) and (90, 120).  
      If \code{y} is simple (time, status) survival observation
      intervals implicitly start at 0, so in that case any duplicate
      identifiers will generate an overlap.}
   \item{gap}{one or more gaps in a subject's timeline; where they are
      in the same state at their return as when they left.}
    \item{jump}{a hole in a subject's timeline, where they are in one
      state at the end of the prior interval, but a new state in the
      at the start subsequent interval.}
    \item{teleport}{two adjacent intervals for a subject, with the
      first interval ending in one state and the subsequent interval
      starting in another.  They have instantaneously changed states
      in 0 units of time.}
    \item{duplicate}{not currently used}
     }
  The total number of occurences of each is present in the \code{flags}
  vector. Optional components give the location and identifiers of the
  flagged observations.
  The \code{Surv} function has already flagged any 0 length intervals as errors.
}
\value{
  a list with components
  \item{states}{the vector of possible states, a union of what appears
  in the Surv object and \code{istate}, with initial states first}
  \item{transitions}{a matrix giving the count of transitions from one
    state to another}
  \item{statecount}{table of the number of visits per state, e.g., 18
    subjects had 2 visits to the "infection" state}
  \item{flags}{a vector giving the counts of each check}
  \item{istate}{a constructed istate that best satisfies all the checks}
  \item{overlap}{a list with the row number and id of overlaps (not
    present if there are no overlaps)}
  \item{gaps}{a list with the row number and id of gaps (not present if
    there are no gaps)}
  \item{teleport}{a list with the row number and id of inconsistent
    rows (not present if there are none)}
  \item{jumps}{a list with the row number and id of jumps (not present
    if there are no jumps)}
}
\note{
For data sets with time-dependent covariates, a given subject will often
have intermediate rows with a status of `no event at this time', coded
as the first level of the factor variable in the Surv() call.
For instance a subject who started in state 'a' at time 0, transitioned to state
'b' at time 10, had a covariate \code{x} change from 135 to 156 at time
20, and a final transition to state 'c' at time 30.
The response would be \code{Surv(c(0, 10, a), c(10, 20, censor),
  c(20,0,c))} where the state variable is a factor with levels of censor,
  a, b, c.
The state variable records \emph{changes} in state, and there was no
change at time 20.
The \code{istate} variable would be (a, b, b); it contains the \emph{current}
state, and the value is unchanged when status = censored.
(It behaves like a \code{tdc} variable from \code{tmerge}).

The intermediate time above is not actually censoring, i.e., a point at
which follow-up for the observation ceases.
The 'censor' label is traditional, but 'none' may be a more accurate choice.

When there are intermediate observations \code{istate} is not
simply a lagged version of the state, and may be more challenging to
create.
One approach is to let \code{survcheck} do the work: call it with
an \code{istate} argument that is correct for the first row of each
subject, or no \code{istate} argument at all, and then insert the
returned value into a data frame.
}
\keyword{ survival }

