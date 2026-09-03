\name{gastric}
\alias{gastric}
\docType{data}
\title{
Gastric cancer data set
}
\description{
Gastric cancer data used the the Moreau paper
}
\usage{data("gastric")}
\format{
  A data frame with 90 observations on the following 3 variables.
  \describe{
    \item{\code{time}}{follow-up time}
    \item{\code{status}}{0=censored, 1= event}
    \item{\code{treat}}{a factor with levels \code{chemo} and
      \code{chemo+radiation}}
  }
}
\note{
  The usege above assumes that the survival package is attached via
  \code{library(survival)} or via another package that depends on
  survival, e.g., mstate.  If not, then either \code{survival::gastric}
  or \code{data(cancer, package="survival") can be used}.
}
\references{
Moreau,T., O'Quigley,J., and Mesbah M(1985). A global goodness-of-fit
statistic for the  proportional hazards model,
         Appl. Statist., 34, 212-218 

}
\examples{
gfit <- coxph(Surv(time, status) ~ treat, gastric)
cox.zph(gfit)
}
\keyword{datasets}
