\name{survfit.object}
\alias{survfit.object}
\title{
  Survival Curve Object 
}
\description{
This class of objects is returned by the \code{survfit} class of functions
to represent a fitted survival curve. 

Objects of this class have methods for the functions \code{print},
\code{summary}, \code{plot}, \code{points} and \code{lines}. The
\code{\link{print.survfit}} method does more computation than is typical
for a print method and is documented on a separate page.
Class of objects that represent a fitted survival curve. 
}
\section{Structure}{
The following components must be included in a legitimate 
\code{survfit} 
object. 
}
\arguments{
\item{n}{
total number of subjects in each curve.
}
\item{time}{
the time points at which the curve has a step. 
}
\item{n.risk}{
the number of subjects at risk at t. 
}
\item{n.event}{
the number of events that occur at time t. 
}
\item{n.enter}{
for counting process data only,
the number of subjects that enter at time t. 
}
\item{n.censor}{
for counting process data only,
the number of subjects who exit the risk set,
without an event,  at time t. 
(For right censored data, this number can be computed from the successive
values of the number at risk).
}
\item{surv}{
the estimate of survival at time t+0. 
This may be a vector or a matrix. 
}
\item{std.err}{
the standard error of the cumulative hazard or -log(survival). 
}
\item{upper}{
upper confidence limit for the survival curve. 
}
\item{lower}{
lower confidence limit for the survival curve. 
}
\item{strata}{
if there are multiple curves, this component gives the number of elements of 
the \code{time} etc. vectors corresponding to the first curve,
the second curve, 
and so on.  The names of the elements are labels for the curves. 
}
\item{start.time}{
the value specified for the \code{start.time} argument,
if it was used in the call.  
}
\item{n.all}{
for counting process data, and any time that the 
\code{start.time} argument was used, this contains
the total number of observations that were available.
Not all may have been used in creating the curve, in which case this value
will be larger than \code{n} above.
of observations that were available
}
\item{conf.type}{
the approximation used to compute the confidence limits. 
}
\item{conf.int}{
the level of the confidence limits, e.g. 90 or 95\%. 
}
\item{na.action}{
the returned value from the na.action function, if any.  It will be used 
in the printout of the curve, e.g., the number of observations deleted due 
to missing values. 
}
\item{call}{
an image of the call that produced the object. 
}
\item{type}{
type of survival censoring. 
}
}
\section{Subscripts}{
Survfit objects that contain multiple survival curves can be subscripted. 
This is most often used to plot a subset of the curves. 
Usually a single subscript will be used.  In one particular case --
survival curves for multiple covariate values, from a Cox model that includes 
a \code{strata} statement -- there is a matrix of 
curves and 2 subscripts may be used. 
(In this case \code{summary.survfit} will 
also print the data as a matrix). 
}
\seealso{
\code{\link{plot.survfit}},  
\code{\link{summary.survfit}},
\code{\link{print.survfit}},
\code{\link{survfit}}.   
}
\keyword{survival}
