\name{Surv}
\alias{Surv}
\alias{is.Surv}
\alias{[.Surv}
\title{
  Create a Survival Object 
}
\description{
Create a survival object, usually used as a response variable in a model 
formula. Argument matching is special for this function, see Details below.
}
\usage{
Surv(time, time2, event,
    type=c('right', 'left', 'interval', 'counting', 'interval2'),
    origin=0)
is.Surv(x)
}
\arguments{
  \item{time}{
    for right censored data, this is the follow up time.  For interval
    data, the first argument is the starting time for the interval. 
  }
  \item{event}{
    The status indicator, normally 0=alive, 1=dead.  Other choices are
    \code{TRUE}/\code{FALSE} (\code{TRUE} = death), 1/2 (2=death), or
    a factor variable.
    For interval censored data, the status indicator is 0=right censoed,
    1=event at \code{time}, 2=left censored, 3=interval censored.
    For multiple endpoint data the event variable will always be a factor,
    whose first level is treated as censoring, or more formally
    "no transtition at this time point".
    There is no constraint on the labels of the factor.
    Although unusual, the event indicator can be omitted, in which case
    all subjects are assumed to have an event.
  }
  \item{time2}{
    ending time of the interval for interval censored  or counting
    process data only.  Intervals are assumed to be open on the left and
    closed on the right, \code{(start, end]}.  For counting process
    data, \code{event} indicates whether a transtion to another state
    occurred at the end of the interval.
  }
  \item{type}{
    character string specifying the type of censoring. Possible values
    are \code{"right"}, \code{"left"}, \code{"counting"},
    \code{"interval"}, \code{"interval2"}.  The default is multi-state
    if \code{event} is a factor, counting process if \code{time2} is
    present, or right censored, in that order. 
  }
  \item{origin}{
    for counting process data, the hazard function origin.  This option
    was intended to be used in conjunction with a model containing
    time dependent
    strata in order to align the subjects properly when they cross over
    from one strata to another, but it has rarely proven useful and is
    depricated.}
  \item{x}{
    any R object.
  }
}

\value{
  An object of class \code{Surv}.  There are methods for \code{print},
  \code{is.na}, and subscripting survival objects.   \code{Surv} objects
  are implemented as a matrix of 2 or 3 columns that has further
  attributes. These include the type (left censored, right censored,
  counting process, etc.) and labels for the states for multi-state
  objects.  Any attributes of the input arguments are also preserved
  in \code{inputAttributes}.  This may be useful for other packages that
  have attached further information to data items such as labels; none
  of the routines in the survival package make use of these
  values, however.

  In the case of \code{is.Surv}, a logical value \code{TRUE} if \code{x}
  inherits from class \code{"Surv"}, otherwise an \code{FALSE}.
}
\details{
   Interval censored data can be represented in two ways.  For the first
   use \code{type = "interval"} and the codes shown above.  In that usage the
   value of the \code{time2} argument is ignored unless event=3.
   The second approach is to think of each observation as a time
   interval with (-infinity, t2) for left censored, (t1, infinity) for
   right censored, (t,t) for exact and (t1, t2) for an interval.
   This is the approach used for type = interval2.  Infinite values can
   be represented either by actual infinity (Inf) or NA.
   The second form has proven to be the more useful one.

Presently, the only methods allowing interval censored data are the 
parametric models computed by \code{survreg} and survival curves
computed by \code{survfit}; for both of these, 
the distinction between open and closed intervals
is unimportant.  
The distinction is important for counting process data and 
the Cox model. 

The function tries to distinguish between the use of 0/1 and 1/2 coding for 
censored data via the condition 
\code{if (max(status)==2)}. 
If 1/2 coding is used and all the subjects are censored, it will 
guess wrong.
In any questionable case it is safer to use logical coding,
e.g., \code{Surv(time, status==3)} would indicate that '3' is
the code for an event.
For multi-state survival the status variable will be a factor, whose
first level is assumed to correspond to censoring.
If there are no uncensored subjects, ensure that the levels of the
factor still contain the censored level!

Surv objects can be subscripted either as a vector, e.g.
\code{x[1:3]} using a single subscript,
in which case the \code{drop} argument is ignored and the result will be
a survival object; 
or as a matrix by using two subscripts.
If the second subscript is missing and \code{drop=F}
(the default),
the result of the subscripting will be a Surv object, e.g., 
\code{x[1:3,,drop=F]},
otherwise the result will be a matrix (or vector), in accordance with
the default behavior for subscripting matrices. 
}

\note{
  The use of 1/2 coding for status is an interesting historical artifact.
  For data contained on punch cards, IBM 360 Fortran treated blank as a zero,
  which led to a policy within the Mayo Clinic section of Biostatistics to never
  use "0" as a data value, since one could not distinguish it from a
  missing value (blank).
  Policy became habit, as is often the case, and the use of 1/2 coding for
  no/yes variables, including death, endured long after the demise of
  the punch cards that had sired the practice.
  At the time \code{Surv} was written many Mayo data sets still used this
  obsolete convention, e.g., the \code{lung} data set found in the package.
}   
  
\seealso{
\code{\link{coxph}},  
\code{\link{survfit}},  
\code{\link{survreg}}, \code{\link{lung}}.   
}
\examples{
with(aml, Surv(time, status))
survfit(Surv(time, status) ~ ph.ecog, data=lung)
Surv(heart$start, heart$stop, heart$event) 
}
\keyword{survival}

