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
  An object of class \code{summary.coxph}, with components:
  \item{n, nevent}{number of observations and number of events,
    respectively, in the fit}
  \item{loglik}{the log partial likelihood at the initial and final
    values}
  \item{coefficients}{a matrix with one row for each coefficient, and
    columns containing the coefficient, the hazard ratio exp(coef),
    standard error, Wald statistic, and P value.}
  \item{conf.int}{a matrix with one row for each coefficient, containing
    the confidence limits for exp(coef)}
  \item{logtest, sctest, waldtest}{the overall likelihood ratio, score,
    and Wald test statistics for the model}
  \item{concordance}{the concordance statistic and its standard error}
  \item{used.robust}{whether an asymptotic or robust variance was used}
  \item{rsq}{an approximate R^2 based on Nagelkirke (Biometrika 1991).}
  \item{fail}{a message, if the underlying coxph call failed}
  \item{call}{a copy of the call}
  \item{na.action}{information on missing values}
}
\note{
  The pseudo r-squared of Nagelkirke is attractive because it is simple,
  but further work has shown that it has poor properties and it is now
  deprecated.  The value is no longer printed by default, and will
  eventually be removed from the object.
}
\seealso{
  \code{\link{coxph}}, \code{\link{print.coxph}}  
}
\examples{
fit <- coxph(Surv(time, status) ~ age + sex, lung) 
summary(fit)
}
\keyword{survival}
