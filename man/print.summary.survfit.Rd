\name{print.summary.survfit}
\alias{print.summary.survfit}
\title{
  Print Survfit Summary 
}
\description{
Prints the result of \code{summary.survfit}. 
}
\usage{
\method{print}{summary.survfit}(x, digits = max(options() $digits-4, 3), ...)
}
\arguments{
\item{x}{
an object of class \code{"summary.survfit"}, which is the result of the 
\code{summary.survfit} function. 
}
\item{digits}{
the number of digits to use in printing the numbers. 
}
  \item{\dots}{for future methods}
}
\value{
\code{x}, with the invisible flag set to prevent printing. 
}
\section{Side Effects}{
prints the summary created by \code{summary.survfit}. 
}
\seealso{
\code{\link{options}},  \code{\link{print}},  \code{\link{summary.survfit}}.   
}
\keyword{print}
% docclass is function
% Converted by Sd2Rd version 37351.
