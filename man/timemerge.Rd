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
tmerge(data, id, start = "start", end, status, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{the primary data set, to which new variables and/or
    observation will be added}

  \item{id}{subject identifier.  For an initial call it identifies the
    variable name of the id variable in the main data set.}
  \item{start}{the name of the variable for the start of time intervals.
    If missing a new variable \code{start} will be added, initialized to
    zero.}

  \item{end}{the variable name that identifies the end of each time interval}

  \item{status}{the variable name that describes status.  It must have
    values of either a TRUE/FALSE or 0/1}

  \item{\dots}{the variables to be added}
}
\details{
 The program is usually run in multiple passes, the first of which
 defines the basic structure, and subsequent ones that add new
 variables to that structure.  The first call will use the
 \code{data}, \code{id}, \code{start},\code{end} and \code{status}
 arguments, while subsequent ones use the \code{data}, \code{id} and
 optional arguments.

 There are 5 types of optional arguments.
 \itemize{
   \item{newname = count(y)}{A new variable will created which contains
   an event count.  The argument \code{y} is assumed to be on the
   scale of the start and end time, and each instance decribes the
   occurent of a "condition" at that time.  A given subject can have
   multiple events at the same time, which will lead to counts of
   other than 0 or 1.}
 \item{newname = cumcount(y)}{Similar to count, except that the event
   count is accumulated over time for each subject.  If x were the first
   date of diagnosis for diabetes, say, then a count variable would
   denote the interval in which diagnosis occured where cumcount
   will denote a current or prior diagnosis.}
 \item{newname = event(y)}{Mark that event occured at time y.}
 \item{newname = cumevent(y)}{Cumulative events}.
 \item{newname = variable(y, x)}{A variable changes its value at time
   y to the new value x}.
 }

 Say that a subject had an interval of observation from age 17 to 38,
 denoted as (17, 38] and that a count occurs at age 24.  A count is a
 predictor which is assumed from the time it occured to the end of
 the interval.  The updated data set will have intervals of (17,24] and
 (24, 38] with a count of 0 for the first and 1 for the second,
 assuming no other occurences for this subject at that age.
 The \code{cumcount} and \code{variable} markers act similarly.
 An event is an outcome, so if coded as an event the said occurence
 would be placed in the (17,24] interval marking that this interval
 has an event at its endpoint.
}
\value{a data frame with two extra attributes \code{tname} and
  \code{tcount}.
  The first contains the names of the key variables; it's persistence
  from call to call allows the user to avoid constantly retyping them.
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
  A \code{count} or \code{variable} code wil count into
  this interval, but not an \code{event} code.  Symmetrically an \code{event}
  occurence at 35 would count in the (17,35] interval, but not one of type
  \code{count}.  The last case is where the main data set has touching
  intervals for a subject, e.g. (17, 28] and (28,35] and a new occurence
  lands at the join.  Events will go to the earlier interval and counts
  to the latter one.

  It is always wise to look at \code{attr(data, 'tcount')} after each
  step of a data set build to avoid surprises.
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
newdata <- with(jasa, data.frame(id= 1:nrow(jasa), 
                  age=     as.numeric(accept.dt - birth.dt)/365.25,
                  futime = as.numeric(fu.date - accept.dt),
                  fustat = fustat, surgery=surgery, mismatch=mismatch))
newdata$futime <- pmax(.5, newdata$futime) # Subject who died on day 0
txtime <- with(jasa, as.numeric(tx.date - accept.dt))
txtime[38] <- txtime[38] - .1

newdata <- tmerge(newdata, id="id", end="futime", status="fustat")
newdata <- tmerge(newdata, id= 1:nrow(jasa), tx = count(txtime))

fit <- coxph(Surv(start, futime, fustat) ~ tx + age, data=newdata)

attr(newdata, "tcount")
# Two subjects transplanted on the day of entry, the "front edge" of
#  their follow-up interval
}
\keyword{ survival }

