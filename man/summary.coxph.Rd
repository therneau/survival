\name{summary.coxph}
\alias{summary.coxph}
\title{
  Summary method for Cox models 
}
\description{
Produces a summary of a fitted coxph model 
}
\usage{
\method{summary}{coxph}(object, conf.int=0.95, scale=1,...)
}
\arguments{
  \item{object}{
    the result of a coxph fit 
  }
  \item{conf.int}{
    level for computation of the confidence intervals. 
    If set to FALSE no confidence intervals are printed 
  }
  \item{scale}{
    vector of scale factors for the coefficients, defaults to 1. 
    The printed coefficients, se, and confidence intervals will be
    associated with one scale unit. 
  }
  \item{\dots}{for future methods}
}
\value{
An object of class \code{summary.coxph}. 
}
\seealso{
  coxph, print.coxph  
}
\examples{
fit <- coxph(Surv(time, status) ~ age + sex, lung) 
summary(fit)
}
\keyword{survival}
