\name{summary.survfit}
\alias{summary.survfit}
\title{
  Summary of a Survival Curve 
}
\description{
  Returns a list containing the survival curve, confidence limits for the 
  curve, and other information. 
}

\usage{
  \method{summary}{survfit}(object, times=, censored=FALSE, scale=1,
  extend=FALSE, rmean=getOption('survfit.rmean'), ...)
  }
\arguments{
    \item{object}{
      the result of a call to the \code{survfit} function. 
    }
    \item{times}{
    vector of times; 
    the returned matrix will contain 1 row for each time. 
    This must be in increasing order and missing values are not allowed. 
    If \code{censored=T}, the default \code{times} vector contains all
    the unique times in \code{fit}, otherwise
    the default \code{times} vector uses only the event (death) times. 
    }
    \item{censored}{
    logical value:  should the censoring times be included in the output?
    This is ignored if the \code{times} argument is present. 
    }
    \item{scale}{
    numeric value to rescale the survival time, e.g., if the input data to 
    \code{survfit} were in 
    days, \code{scale = 365.25} would scale the output to years. 
    }
    \item{extend}{
    logical value: if TRUE, prints information for all specified \code{times}, 
    even if there are no subjects left at the end of the specified
    \code{times}. 
    This is only valid if the \code{times} argument is present. 
  }
  \item{rmean}{Show restricted mean: see
  \code{\link{print.survfit}} for details}
  \item{\dots}{for future methods}
  }

\value{
a list with the following components: 

\item{surv}{
the estimate of survival at time t+0. 
}
\item{time}{
the timepoints on the curve. 
}
\item{n.risk}{
the number of subjects at risk at time t-0 
(but see the comments on weights in the \code{survfit} help file). 
}
\item{n.event}{
if the \code{times} argument is missing, then this column is the number of 
events that occurred at time t. 
Otherwise, it is the cumulative number of events that have occurred 
since the last time listed until time t+0. 
}
\item{n.entered}{
This is present only for counting process survival data.
If the \code{times} argument is 
missing, this column is the number of subjects that entered at time t. 
Otherwise, it is the cumulative number of subjects that have entered 
since the last time listed until time t.  
}
\item{n.exit.censored}{

if the \code{times} argument is 
missing, this column is the number of subjects that left without an 
event at time t. 
Otherwise, it is the cumulative number of subjects that have left 
without an event 
since the last time listed until time t+0.  
This is only present for counting process survival data.
}
\item{std.err}{
the standard error of the survival value. 
}
\item{conf.int}{
level of confidence for the confidence intervals of survival. 
}
\item{lower}{
lower confidence limits for the curve. 
}
\item{upper}{
upper confidence limits for the curve. 
}
\item{strata}{
indicates stratification of curve estimation.  
If \code{strata} is not \code{NULL}, 
there are multiple curves in the result and the \code{surv}, \code{time}, \code{n.risk}, etc.  
vectors will contain multiple curves, pasted end to end.  
The levels of \code{strata} (a factor) are the labels for the curves. 
}
\item{call}{
the statement used to create the \code{fit} object. 
}
\item{na.action}{
same as for \code{fit}, if present. 
}
\item{table}{
table of information that is returned from \code{print.survfit} function. 
}
\item{type}{
type of data censoring.  Passed through from the fit object. 
}
}
\seealso{
\code{\link{survfit}},  \code{\link{print.summary.survfit}}   
}
\examples{
summary( survfit( Surv(futime, fustat)~1, data=ovarian))
summary( survfit( Surv(futime, fustat)~rx, data=ovarian))
}
\keyword{survival}
