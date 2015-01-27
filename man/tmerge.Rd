\name{tmerge}
\alias{tmerge}
\title{Time based merge for survival data}
\description{
  A common task in survival analysis is the creation of start,stop data
  sets which have multiple intervals for each subject, along with the
  covariate values that apply over that interval.  This function aids
  in the creation of such data sets.
}
\usage{
tmerge(data, id,\dots, tname)
}
\arguments{
  \item{data}{the primary data set, to which new variables and/or
    observation will be added}

  \item{id}{subject identifier}

  \item{\dots}{operations that add new variables or intervals, see below}
  \item{tname}{optional, a character vector or list that identifies the
    key variable names in \code{data}.  These are the variables for the
  start and end of the time intervals, the subject identifier, and
  optionally a status variable.  The default name for the start time of
 the intervals is "tstart" if none is given.}
}
\details{
 The program is usually run in multiple passes, the first of which
 defines the basic structure, and subsequent ones that add new
 variables to that structure.  The first call will use the
 \code{data}, \code{id}, \code{start},\code{end} and \code{status}
 arguments, while subsequent ones use the \code{data}, \code{id} and
 optional arguments.  For a more complete explanation of how this
 routine works refer to the vignette on time-dependent variables.

 There are 4 types of optional arguments: a count, cumulative count,
 event, or cumulative event.  Counts are time dependent covariates
 which occur before an event, events are outcomes
 \itemize{
   \item{newname = count(y, x)}{A new variable will created which contains
   an event count.  The argument \code{y} is assumed to be on the
   scale of the start and end time, and each instance decribes the
   occurent of a "condition" at that time.
   The second argument \code{x} is optional. In the usual case where
   \code{x} is missing the count variable starts at 0 for each subject
   and becomes 1 at the time of the event; if \code{x} is present the
   count is set to the value of \code{x}.
   If a given subject has multiple rows of data with the same time value
   the sum of those rows will be assigned (number of rows if \code{x} is
   missing).  

 \item{newname = cumcount(y,x)}{Similar to count, except that the event
   count is accumulated over time for each subject.}

 \item{newname = event(y,x)}{Mark that event occured at time y.
   In the ususal case that \code{x} is missing, the new 0/1 variable
   will be similar to the 0/1 status variable of a survival time, and
   that is in fact how it will normally be used.  For multiple types
   of endpoints the \code{x} variable can be used encode the type of
   event.
 }
 \item{newname = cumevent(y,x)}{Cumulative events}.
 }
}
 Say that a subject had an interval of observation from age 17 to 38,
 denoted as (17, 38] and that a count occurs at age 24.  A count is a
 predictor which is assumed to apply from the time it occured to the end of
 the interval.  The updated data set will have intervals of (17,24] and
 (24, 38] with a count of 0 for the first interval and 1 for the second,
 assuming no other occurences for this subject at exactly time 24.
 An event is an outcome, so if coded as an event the said occurence
 would be placed in the (17,24] interval, with the new variable
 marking that this interval finished with an event.
}
\value{a data frame with two extra attributes \code{tname} and
  \code{tcount}.
  The first contains the names of the key variables; it's persistence
  from call to call allows the user to avoid constantly reentering the
  \code{tname} argument.
  The tcount variable contains counts of the match types.
  New time values that occur before the first interval for a subject
  are "early", those after the last interval for a subject are "late",
  and those that fall into a gap are of type "gap".
  For type \code{count} and \code{cumcount} early events are counted in
  the first interval for a subject, to make it easier to capture the
  common "any prior history of" type of variable.  Otherwise early,
  late, and gap observations are ignored, as are those with a missing
  value for the occurence time.

  The most common type will usually be "mid", for those new times that
  fall inside an existing interval and cause it to be split into two.
  Observations that fall exactly on the edge of an interval are counted
  as "front", "back" or "tie".  The first corresponds for instance
  to an occurence at 17 for someone with an interval (17, 35] who is not
  at risk just before time 17.
  A \code{count}  will be added to this interval
  but not an \code{event}.  Symmetrically an \code{event}
  occurence at 35 would count in the (17,35] interval, but a \code{count}
  would not.  The last case is where the main data set has touching
  intervals for a subject, e.g. (17, 28] and (28,35] and a new occurence
  lands at the join.  Events will go to the earlier interval and counts
  to the latter one.

  It is often wise to look at \code{attr(data, 'tcount')} after each
  step of a data set build to avoid surprises.

  These extra attributes are ephemeral, and will be discarded
  if the dataframe is modified in any way.  This is intentional.
}
\author{Terry Therneau}
\seealso{\code{\link{neardate}}}

\examples{
# The data set jasa contains the famous Stanford Heart Transplant data
#  set, as it appeared in Crowley and Hu, JASA 72:27-36, 1971.
# Two special cases need to be dealt with:
#  subject 15 died on day 0 which leads to an illegal (0,0] interval,
#     make them die on day 0.5 instead
#  subject 38 dies on the day of transplant, make tx happen "earlier in
#     the day" (before death) by subtracting .1 from their transplant day
#
temp <- with(jasa, data.frame(id= 1:nrow(jasa), 
             age=     as.numeric(accept.dt - birth.dt)/365.25,
             futime = as.numeric(fu.date - accept.dt),
             fustat = fustat, surgery=surgery, mismatch=mismatch))
temp$futime <- pmax(.5, temp$futime) # Subject who died on day 0
txtime <- with(jasa, as.numeric(tx.date - accept.dt))
txtime[38] <- txtime[38] - .1

# Now create the start/stop data set.  First call is setup, second the work
jasa2 <- tmerge(temp, tname=c(id="id", end="futime", status="fustat"))
jasa2 <- tmerge(jasa3, id= 1:nrow(jasa), transplant = count(txtime))

fit <- coxph(Surv(tstart, futime, fustat) ~ transplant + age, data=jasa2)

attr(jasa2, "tcount")
# Shows two subjects transplanted on the day of entry, the "front edge" of
#  their follow-up interval
}
\keyword{ survival }

