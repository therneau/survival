\name{survexp.object}
\alias{survexp.object}
\title{
  Expected Survival Curve Object 
}
\description{
This class of objects is returned by the \code{survexp} class of functions
to represent a fitted survival curve. 

Objects of this class have methods for \code{summary}, and inherit
the \code{print}, \code{plot}, \code{points} and \code{lines}
methods from \code{survfit}.
}
\section{Structure}{
The following components must be included in a legitimate 
\code{survfit} 
object. 
}
\arguments{
\item{surv}{
the estimate of survival at time t+0. 
This may be a vector or a matrix. 
}
\item{n.risk}{
the number of subjects who contribute at this time.
}
\item{time}{
the time points at which the curve has a step. 
}
\item{std.err}{
the standard error of the cumulative hazard or -log(survival). 
}
\item{strata}{
if there are multiple curves, this component gives the number of elements of 
the \code{time} etc. vectors corresponding to the first curve,
the second curve, 
and so on.  The names of the elements are labels for the curves. 
}
\item{method}{the estimation method used.  One of "Ederer", "Hakulinen",
  or "conditional".}
\item{na.action}{
the returned value from the na.action function, if any.  It will be used 
in the printout of the curve, e.g., the number of observations deleted due 
to missing values. 
}
\item{call}{
an image of the call that produced the object. 
}
}
\section{Subscripts}{
Survexp objects that contain multiple survival curves can be subscripted. 
This is most often used to plot a subset of the curves. 
}
\section{Details}{
In expected survival each subject from the data set is matched to a
hypothetical person from the parent population, matched on the
characteristics of the parent population.
The number at risk printed here is the number of those hypothetical
subject who are still part of the calculation.
In particuarl, for the Ederer method all hypotheticals are retained for
all time, so \code{n.risk} will be a constant.
}
\seealso{
\code{\link{plot.survfit}},  
\code{\link{summary.survexp}},
\code{\link{print.survfit}},
\code{\link{survexp}}.   
}
\keyword{survival}
