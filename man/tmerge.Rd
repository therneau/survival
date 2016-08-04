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
  \item{data2}{second data set in which the other arguments
    will be found}
  \item{id}{subject identifier}

  \item{\dots}{operations that add new variables or intervals, see
    below}
  \item{tstart}{optional variable to define the valid time range for
    each subject, only used on an initial call}
   \item{tstop}{optional variable to define the valid time range for
    each subject, only used on an initial call}
  
  \item{options}{a list of options.  Valid ones are idname, tstartname,
    tstopname, delay and na.rm.  See the explanation below.}
}
\details{
 The program is often run in multiple passes, the first of which
 defines the basic structure, and subsequent ones that add new
 variables to that structure.  For a more complete explanation of how this
 routine works refer to the vignette on time-dependent variables.

 There are 4 types of optional arguments: a time dependent covariate
 (tdc), cumulative count (cumtdc), event (event) or cumulative event
 (cumevent).
 Time dependent covariates change their values before an event,
 events are outcomes.
 \itemize{
   \item{newname = tdc(y, x)}{ A new time dependent covariate
     variable will created. 
   The argument \code{y} is assumed to be on the
   scale of the start and end time, and each instance decribes the
   occurent of a "condition" at that time.
   The second argument \code{x} is optional. In the  case where
   \code{x} is missing the count variable starts at 0 for each subject
   and becomes 1 at the time of the event; if \code{x} is present the
   count is set to the value of \code{x}.
   If the option \code{na.rm=TRUE} missing values of \code{x} are
   first removed.

 \item{newname = cumtdc(y,x)}{ Similar to tdc, except that the event
   count is accumulated over time for each subject.}

 \item{newname = event(y,x)}{ Mark an event at time y.
   In the ususal case that \code{x} is missing the new 0/1 variable
   will be similar to the 0/1 status variable of a survival time.
 }
 \item{newname = cumevent(y,x)}{ Cumulative events}.
 }
}

The function adds three new variables to the output data set: 
\code{tstart}, \code{tstop}, and \code{id}.  The \code{options} argument
can be used to change these names.
The \code{na.rm} option affects creation of time-dependent covariates.
Should a data row in \code{data2} that has a missing value for the
variable be ignored (na.rm=FALSE, default) or should it generate an
observation with a value of NA?  The default value leads to
"last value carried forward" behavior.
The \code{delay} option causes a time-dependent covariate's new
value to be delayed, see the vignette for an example.
}

\value{a data frame with two extra attributes \code{tname} and
  \code{tcount}.
  The first contains the names of the key variables; it's persistence
  from call to call allows the user to avoid constantly reentering the
  \code{options} argument.
  The tcount variable contains counts of the match types.
  New time values that occur before the first interval for a subject
  are "early", those after the last interval for a subject are "late",
  and those that fall into a gap are of type "gap".
  All these are are considered to be outside the specified time frame for the
  given subject.  An event of this type will be discarded.  A
  time-dependent covariate value will be applied to later intervals but
  will not generate a new time point in the output.

  The most common type will usually be "within", for those new times that
  fall inside an existing interval and cause it to be split into two.
  Observations that fall exactly on the edge of an interval but within the
  (min, max] time for a subject are counted
  as being on a "leading" edge, "trailing" edge or "boundary".
  The first corresponds for instance
  to an occurence at 17 for someone with an intervals of (0,15] and (17, 35].
  A \code{tdc} at time 17  will affect this interval
  but an \code{event} at 17 would be ignored.  An \code{event}
  occurence at 15 would count in the (0,15] interval.
  The last case is where the main data set has touching
  intervals for a subject, e.g. (17, 28] and (28,35] and a new occurence
  lands at the join.  Events will go to the earlier interval and counts
  to the latter one.  A last column shows the number of additions where
  where the id and time point were identical.

  These extra attributes are ephemeral and will be discarded
  if the dataframe is modified in any way.  This is intentional.
}
\author{Terry Therneau}
\seealso{\code{\link{neardate}}}

\examples{
# The pbc data set contains baseline data and follow-up status
# for a set of subjects with primary biliary cirrhosis, while the
# pbcseq data set contains repeated laboratory values for those
# subjects.  
# The first data set contains data on 312 subjects in a clinical trial plus
# 106 that agreed to be followed off protocol, the second data set has data
# only on the trial subjects.
temp <- subset(pbc, id <= 312, select=c(id:sex, stage)) # baseline data
pbc2 <- tmerge(temp, temp, id=id, endpt = event(time, status))
pbc2 <- tmerge(pbc2, pbcseq, id=id, ascites = tdc(day, ascites),
               bili = tdc(day, bili), albumin = tdc(day, albumin),
               protime = tdc(day, protime), alk.phos = tdc(day, alk.phos))

fit <- coxph(Surv(tstart, tstop, endpt==2) ~ protime + log(bili), data=pbc2)
}
\keyword{ survival }

