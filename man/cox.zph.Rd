\name{cox.zph}
\alias{cox.zph}
\alias{[.cox.zph}
\alias{print.cox.zph}
\title{
  Test the Proportional Hazards Assumption of a Cox Regression
}
\description{
Test the proportional hazards assumption for a Cox regression model fit 
(\code{coxph}). 
}
\usage{
cox.zph(fit, transform="km", global=TRUE)
}
\arguments{
  \item{fit}{
    the result of fitting a Cox regression model, using the \code{coxph} function. 
  }
  \item{transform}{
    a character string specifying how the survival times should be transformed 
    before the test is performed. 
    Possible values are \code{"km"}, \code{"rank"}, \code{"identity"} or  a 
    function of one argument. 
  }
  \item{global}{
    should a global chi-square test be done, in addition to the 
    per-variable tests. 
  }
}
\value{
an object of class \code{"cox.zph"}, with components: 

\item{table}{
a matrix with one row for each variable, and optionally a last row for 
the global test. 
Columns of the matrix contain the correlation coefficient between transformed 
survival time and the scaled Schoenfeld residuals, a chi-square, 
and the two-sided p-value. 
For the global test there is no appropriate correlation, so an NA is 
entered into the matrix as a placeholder. 
}
\item{x}{
the transformed time axis. 
}
\item{y}{
the matrix of scaled Schoenfeld residuals.  There will be one column per 
variable and one row per event.  The row labels contain the original event 
times (for the identity transform, these will be the same as \code{x}).
}
\item{call}{
the calling sequence for the routine. 

The computations require the original \code{x} matrix of the Cox model fit. 
Thus it saves time if the \code{x=TRUE} option is used in \code{coxph}. 
This function would usually be followed by both a plot and a print of 
the 
result. 
The plot gives an estimate of the time-dependent coefficient \code{beta(t)}. 
If the proportional hazards assumption is true, \code{beta(t)} will be a 
horizontal 
line.  The printout gives a test for \code{slope=0}. 
}
}
\references{
P. Grambsch and T. Therneau (1994), 
Proportional hazards tests and diagnostics based on weighted residuals. 
\emph{Biometrika,}
\bold{81}, 515-26. 
}
\seealso{
\code{\link{coxph}},  \code{\link{Surv}}.   
}
\examples{
fit <- coxph(Surv(futime, fustat) ~ age + ecog.ps,  
             data=ovarian) 
temp <- cox.zph(fit) 
print(temp)                  # display the results 
plot(temp)                   # plot curves 
}
\keyword{survival}
