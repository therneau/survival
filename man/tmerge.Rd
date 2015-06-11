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
tmerge(data1, data2,  id,\dots, tstart, tstop, options)
}
\arguments{
  \item{data1}{the primary data set, to which new variables and/or
    observation will be added}
  \item{data2}{optional second data set in which the other arguments
    will be found}
  \item{id}{subject identifier}

  \item{\dots}{operations that add new variables or intervals, see
    below}
  \item{tstart}{optional variable to define the valid time range for
    each subject, only used on an initial call}
   \item{tstop}{optional variable to define the valid time range for
    each subject, only used on an initial call}
  
  \item{options}{a list of options.  Valid ones are id, tstart, and
    tstop, which will be the names of the three mandatory variables
    in the output data.  The other is defer, which sets a numeric
    amount of time before an event when covariate changes are
    disallowed (the are deferred until after the event in that case.)}
}
\details{
 The program is usually run in multiple passes, the first of which
 defines the basic structure, and subsequent ones that add new
 variables to that structure.  For a more complete explanation of how this
 routine works refer to the vignette on time-dependent variables.

 There are 4 types of optional arguments: a time dependent covariate
 (tdc), cumulative count (cumtdc), event (event) or cumulative event
 (cumevent).
 Time dependent covariates change their values before an event,
 events are outcomes.
 \itemize{
   \item{newname = tdc(y, x)}{A new time dependent covariate
     variable will created. 
   The argument \code{y} is assumed to be on the
   scale of the start and end time, and each instance decribes the
   occurent of a "condition" at that time.
   The second argument \code{x} is optional. In the  case where
   \code{x} is missing the count variable starts at 0 for each subject
   and becomes 1 at the time of the event; if \code{x} is present the
   count is set to the value of \code{x}.
   If a given subject has multiple rows of data with the same time value
   the sum of those rows will be assigned.

 \item{newname = cumtdc(y,x)}{Similar to tdc, except that the event
   count is accumulated over time for each subject.}

 \item{newname = event(y,x)}{Mark an event at time y.
   In the ususal case that \code{x} is missing, the new 0/1 variable
   will be similar to the 0/1 status variable of a survival time, and
   that is in fact how it will normally be used.  For multiple types
   of endpoints the \code{x} argument can be used encode the type of
   event.
 }
 \item{newname = cumevent(y,x)}{Cumulative events}.
 }
}
 Say that a subject had an interval of observation from age 17 to 38,
 denoted as (17, 38] and that a marker occurs at age 24.  A tdc
 variable is a
 predictor which is assumed to apply from the time it occured to the end of
 followup for the subject.
 The updated data set will have intervals of (17,24] and
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
  \code{options} arguments.
  The tcount variable contains counts of the match types.
  New time values that occur before the first interval for a subject
  are "early", those after the last interval for a subject are "late",
  and those that fall into a gap are of type "gap".

  The most common type will usually be "within", for those new times that
  fall inside an existing interval and cause it to be split into two.
  Observations that fall exactly on the edge of an interval are counted
  as "leading" edge, "trailing" or "boundary".
  The first corresponds for instance
  to an occurence at 17 for someone with an interval (17, 35] who is not
  at risk just before time 17.
  A \code{tdc} at time 17  will affect this interval
  but not an \code{event}.  Symmetrically an \code{event}
  occurence at 35 would count in the (17,35] interval, but a \code{tdc}
  would not.  The last case is where the main data set has touching
  intervals for a subject, e.g. (17, 28] and (28,35] and a new occurence
  lands at the join.  Events will go to the earlier interval and counts
  to the latter one.

  It is wise to look at \code{attr(data, 'tcount')} after each
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
tdata <- jasa[, -(1:4)]  #leave off the dates, temporary data set
tdata$futime <- pmax(.5, tdata$futime)  # the death on day 0
indx <- with(tdata, which(wait.time == futime))
tdata$wait.time[indx] <- tdata$wait.time[indx] - .5  #the tied transplant
sdata <- tmerge(tdata, tdata, id=1:nrow(tdata), 
                death = event(futime, fustat), 
                trans = tdc(wait.time))
attr(sdata, "tcount")
# Shows two subjects transplanted on the day of entry, the "front edge" of
#  their follow-up interval

fit <- coxph(Surv(tstart, tstop, death) ~ trans + age, data=sdata)
}
\keyword{ survival }

