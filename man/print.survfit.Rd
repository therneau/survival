\name{print.survfit}
\alias{print.survfit}
\title{
  Print a Short Summary of a Survival Curve 
}
\description{
  Print number of observations, number of events, the restricted
  mean survival and its  
standard error, and the median survival with confidence limits for the  
median. 
}
\usage{
\method{print}{survfit}(x, scale=1, digits = max(options()$digits - 4,3),
    print.rmean=getOption("survfit.print.rmean"),
    rmean = getOption('survfit.rmean'),...)
}
\arguments{
\item{x}{
the result of a call to the \code{survfit} function. 
}
\item{scale}{
a numeric value to rescale the survival time, e.g., 
if the input data to survfit were in days, 
\code{scale=365} would scale the printout to years. 
}
\item{digits}{Number of digits to print}
\item{print.rmean,rmean}{Options for computation and display of the
  restricted mean.}
\item{\dots}{for future results}
}
\value{
  x, with the invisible flag set to prevent printing.
  (The default for all print functions in R is to return the object
  passed to them; print.survfit follows the pattern.  If you want to
  capture these printed results for further processing, see the
  \code{table} component of \code{summary.survfit}.)  
}
\section{Side Effects}{
The number of observations, the number of events,  the median survival
with its confidence interval, and optionally the
restricted mean survival (\code{rmean}) and its standard error, are printed.
If there are multiple curves, there is one line of output for each. 
}
\details{
The mean and its variance are based on a truncated estimator.  That is, if the 
last observation(s) is not a death, then the survival curve estimate does not 
go to zero and the mean is undefined.
There are four possible approaches, which are selected by the
\code{rmean} option.
The first is to set the upper limit to a constant, e.g.,\code{rmean=365}.
In this case the reported mean would be the expected number of days, out
of the first 365, that would be experienced by each group.  This is
useful if interest focuses on a fixed period.
Other options are \code{"none"} (no estimate), \code{"common"} and
\code{"individual"}.
For the last of these, the mean is computed as the area under the curve,
over the
range from 0 to the maximum observed time for that curve.
Since the end point is random, values for different curves are not
comparable and the printed standard errors are an underestimate as
they do not take into account this random variation.  This option is
provided mainly for backwards compatability, as this estimate was the
default (only) one in earlier releases of the code.
The \code{"common"} option uses the maximum time for all curves in the
object as a common upper limit for the auc calculation.

The median and its confidence interval are defined by drawing a horizontal 
line at 0.5 on the plot of the survival curve and its confidence bands. 
The intersection of the line with the lower CI band defines the lower limit 
for the median's interval, and similarly for the upper band.  If any of the 
intersections is not a point, then we use the smallest point of intersection, 
e.g., if the survival curve were exactly equal to 0.5 over an interval.
}
\section{References}{
Miller, Rupert G., Jr. (1981).  
\emph{Survival Analysis.}
New York:Wiley, p 71. 
}
\seealso{
\code{\link{summary.survfit}}. 
}
\keyword{survival}
% docclass is function
% Converted by Sd2Rd version 37351.
