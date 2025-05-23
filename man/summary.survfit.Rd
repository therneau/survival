\name{summary.survfit}
\alias{summary.survfit}
\alias{summary.survfitms}
\title{
  Summary of a Survival Curve 
}
\description{
  Returns a list containing the survival curve, confidence limits for the 
  curve, and other information. 
}

\usage{
  \method{summary}{survfit}(object, times, censored=FALSE, scale=1,
  extend=FALSE, rmean=getOption('survfit.rmean'), data.frame=FALSE, dosum, ...)
  \method{summary}{survfitms}(object, times, censored=FALSE, scale=1,
  extend=FALSE, rmean=getOption('survfit.rmean'), data.frame=FALSE, ...)
  }
\arguments{
    \item{object}{
      the result of a call to the \code{survfit} function. 
    }
    \item{times}{
    vector of times; 
    the returned matrix will contain 1 row for each time. 
    The vector will be sorted into increasing order;
    missing values are not allowed. 
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
    This is only used if the \code{times} argument is present. 
  }
  \item{rmean}{Show restricted mean: see
    \code{\link{print.survfit}} for details}
  \item{data.frame}{if TRUE, return the results as a data frame, rather
    than a summary.survfit object}
  \item{dosum}{only applicable if \code{times} is present, see details below}
  \item{\dots}{for future methods}
  }

\value{
  if \code{data.frame = TRUE}, a data frame with columns of time,
  n.risk, n.event, n.censor, surv, cumhaz, strata (if present) and
  data (the row of newdata for survfit.coxph).  Also std.err, std.chaz,
  upper and lower if the curve had se.fit=TRUE.

  if \code{data.frame = FALSE}, a list with the following components: 

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
\section{Details}{
This routine has two uses: printing out a survival curve at specified
time points (often yearly), or extracting the values at specified time
points for further processing.
In the first case we normally want \code{extend=FALSE}, i.e., don't print out
data past the end of the curve.  If the \code{times} option only
contains values beyond the last point in the curve then there is nothing
to print and an error message will result.
For the second usage we often want \code{extend=TRUE}, so that the
results will have a predictable length.
If \code{data.frame = TRUE} then either might be desired.
Be aware, however, that these extended values will often be badly
biased; we are essentialy treating the final censored subjects as
immortal.

The underlying survival object will have a row for each unique event or
censoring time.
When the \code{times} argument contains values not in the data, the
routine can only use a best guess for the number at risk, i.e., the
number at risk at the next event/censoring time.
When the routine is called with counting process data many users are
confused by counts that appear too large.
For example, \code{Surv(c(0,0, 1, 5), c(2, 3, 8, 10), c(1, 0, 1, 0))},
which prints as \code{(0,2] (0, 3+] (1, 8]  (5,10+]}.
Do \code{survfit} followed by \code{summary} with
a request for the values at time 0.
The \code{survfit} object has entries only at times 2, 3, 8, and 10;
there are 3 subjects at risk at time 2, so that is what will be printed
for time 0.

For a printout at fixed times, for example yearly values for a curve,
the printed number of events will by default be the total number of
events that have occured since the prior line of printout, and likewise
for number of censored and number at entry, \code{dosum = TRUE}.
Alternately, the routine can return the number of events/censors/entry at
that time, \code{dosum=FALSE}.
This feature was added at the request of a user who essentially wanted to
use the \code{times} argument as a subscript to pick off selected rows
of the output, e.g., to select survival values corresponding to the
last follow-up times of a new set of observations.
The default for \code{dosum} is TRUE if the \code{times} vector is
strictly increasing and FALSE otherwise.

For a survfitms object replace the \code{surv} component with
\code{pstate}.  Also, a data frame will not include the cumulative
hazard since it has a different multiplicity: one column per transition
rather than one per state.
}
\seealso{
\code{\link{survfit}},  \code{\link{print.summary.survfit}}   
}
\examples{
summary( survfit( Surv(futime, fustat)~1, data=ovarian))
summary( survfit( Surv(futime, fustat)~rx, data=ovarian))
}
\keyword{survival}
