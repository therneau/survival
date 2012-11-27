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
    type=c('right', 'left', 'interval', 'counting', 'interval2', 'mstate'),
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
    \code{"interval"}, \code{"interval2"} or \code{"mstate"}.  
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
  When the \code{type} argument is missing the code assumes a type based
  on the following rules:
  \itemize{
    \item If there are two unnamed arguments, they will match \code{time} and
     \code{event} in that order.  If there are three unnamed arguments
     they match \code{time}, \code{time2} and \code{event}.
     \item If the event variable is a factor then type \code{mstate} is
     assumed.  Otherwise type \code{right} if there is no \code{time2}
     argumement, and type \code{counting} if there is.
   }
   As a consequence the \code{type} argument can usually be omitted.

   When the survival type is "mstate" then the status variable will be
   treated as a factor.  The first level of the factor is taken to
   represent censoring and remaining ones a transition to the given state.
   
    Interval censored data can be represented in two ways.  For the first
   use type = interval and the codes shown above.  In that usage the
   value of the {time2} argument is ignored unless event=3.
   The second approach is to think of each observation as a time
   interval with (-infinity, t) for left censored, (t, infinity) for
   right censored, (t,t) for exact and (t1, t2) for an interval.
   This is the approach used for type = interval2, with NA taking
   the place of infinity.  It has proven to be the more useful.

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

For multi-state survival (type= "mstate") the status variable can have
multiple levels.  The first of these will stand for censoring, and the
others for various event types, e.g., causes of death.  

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

