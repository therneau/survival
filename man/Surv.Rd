\name{Surv}
\alias{Surv}
\alias{is.Surv}
\alias{print.Surv}
\alias{Math.Surv}
\alias{Summary.Surv}
\alias{[.Surv}
\alias{format.Surv}
\alias{as.data.frame.Surv}
\alias{as.character.Surv}
\alias{as.matrix.Surv}
\alias{is.na.Surv}
\alias{Ops.Surv}
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
    \code{TRUE}/\code{FALSE} (\code{TRUE} = death) or 1/2 (2=death). For
    interval censored data, the status indicator is 0=right censored,
    1=event at \code{time}, 2=left censored, 3=interval censored.
    Although unusual, the event indicator can be omitted, in which case
    all subjects are assumed to have an event.
  }
  \item{time2}{
    ending time of the interval for interval censored  or counting
    process data only.  Intervals are assumed to be open on the left and
    closed on the right, \code{(start, end]}.  For counting process
    data, \code{event} indicates whether an event occurred at the end of
    the interval.
  }
  \item{type}{
    character string specifying the type of censoring. Possible values
    are \code{"right"}, \code{"left"}, \code{"counting"},
    \code{"interval"}, or \code{"interval2"}.  The default is
    \code{"right"} or \code{"counting"} depending on whether the
    \code{time2} argument is absent or present, respectively.
  }
  \item{origin}{
    for counting process data, the hazard function origin.  This option
    was intended to be used in conjunction with a model containing
    time dependent
    strata in order to align the subjects properly when they cross over
    from one strata to another, but it has rarely proven useful.}
  \item{x}{
    any R object.
  }
}

\value{
  An object of class \code{Surv}.  There are methods for \code{print},
  \code{is.na}, and subscripting survival objects.   \code{Surv} objects
  are implemented as a matrix of 2 or 3 columns.

  In the case of \code{is.Surv}, a logical value \code{TRUE} if \code{x}
  inherits from class \code{"Surv"}, otherwise an \code{FALSE}.
}
\details{
  Typical usages are
\preformatted{
Surv(time, event)
Surv(time, time2, event, type= )
}
The \code{time},\code{time2} and \code{event} arguments are matched by
position, not by name, so use, eg, \code{Surv(time, dead)} rather than
\code{Surv(time, event=dead)}

In theory it is possible to represent interval censored data without a 
third column containing the explicit status.  Exact, right censored, 
left censored and interval censored observation would be represented as 
intervals of [a,a], (a, infinity), (-infinity,b), and [a,b]  
respectively; each interval is a pair of time points 
within which the event is  
known to have occurred. 

If \code{type="interval2"} then the representation given above is 
assumed, with \code{NA} taking the place of infinity.  
If \code{type="interval"} then \code{event} must be given.  
If \code{event} is \code{0}, \code{1}, or \code{2}, 
the relevant 
information is assumed to be contained in \code{time},
the value in \code{time2} is ignored,
and the second column of the internal representation contains a 
placeholder value.

Presently, the only methods allowing interval censored data are the 
parametric models computed by \code{survreg} and survival curves
computed by \code{survfit}; for both of these, 
the distinction between open and closed intervals
is unimportant.  
The distinction is important for counting process data and 
the Cox model. 

The function tries to distinguish between the use of 0/1 and 1/2 coding for 
left and right censored data using 
\code{if (max(status)==2)}. 
If 1/2 coding is used and all the subjects are censored, it will 
guess wrong.
In any questionable case it is safer to use logical coding,
e.g., \code{Surv(time, status==3)} would indicate that a \code{3} is
the code for an event.

Surv objects can be subscripted either as an object, e.g.
\code{x[1:3]} using a single subscript;
in which case the \code{drop} argument is ignored;
or as a matrix, using two arguments.
If the second subscript is missing and \code{drop=F}
(the default),
the result of the subscripting will be a Surv object, e.g., 
\code{x[1:3,,drop=F]},
otherwise the result will be a matrix (or vector), in accordance with
the default behavior for subscripting matrices. 
}

\seealso{
\code{\link{coxph}},  
\code{\link{survfit}},  
\code{\link{survreg}}.   
}
\examples{
with(lung, Surv(time, status))
Surv(heart$start, heart$stop, heart$event) 
}
\keyword{survival}

