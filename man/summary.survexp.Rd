\name{summary.survexp}
\alias{summary.survexp}
\title{Summary function for a survexp object}
\description{
  Returns a list containing the values of the survival at
  specified times.
}
\usage{
\method{summary}{survexp}(object, times, scale = 1, ...)
}
\arguments{
  \item{object}{
    the result of a call to the \code{survexp} function
}
  \item{times}{
   vector of times; 
    the returned matrix will contain 1 row for each time. 
    Missing values are not allowed.  
}
  \item{scale}{
    numeric value to rescale the survival time, e.g., if the input data to 
    \code{survfit} were in 
    days, \code{scale = 365.25} would scale the output to years. 
}
  \item{\dots}{For future methods}
}
\details{
  A primary use of this function is to retreive survival at fixed time
  points, which will be properly interpolated by the function.
}
\value{
a list with the following components: 

\item{surv}{
the estimate of survival at time t. 
}
\item{time}{
the timepoints on the curve. 
}
\item{n.risk}{
In expected survival each subject from the data set is matched to a
hypothetical person from the parent population, matched on the
characteristics of the parent population.
The number at risk is the number of those hypothetical
subject who are still part of the calculation.
}
}
\author{Terry Therneau}
\seealso{\code{\link{survexp}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
