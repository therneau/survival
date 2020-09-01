\name{royston}
\alias{royston}
\title{Compute Royston's D for a Cox model}
\description{
  Compute the D statistic and R^2 for a coxph model, proposed by
  Royston and Sauerbrei
}
\usage{
royston(fit, newdata, ties = TRUE, adjust = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{a coxph fit}
  \item{newdata}{optional validation data set}
  \item{ties}{make a correction for ties in the risk score}
  \item{adjust}{adjust for possible overfitting}
}
\details{
  The adjustment is based on the ratio
  r= (number of events)/(number of coefficients). For models which have
  sufficient sample size (r>20) the adjustment will be small.
}
\value{a vector containing the value of D, \eqn{R-squared}{R^2_D} and
  the estimated standard error of D.
}
\references{
P. Royston and W. Sauerbrei, A new measure of prognostic separation in
survival data.  Statistics in Medicine 23:723-748, 2004.
}
\examples{
# one of the examples from the paper
pbc2 <- na.omit(pbc)  # no missing values
cfit <- coxph(Surv(time, status==2) ~ age + log(bili) + edema + albumin +
                   stage + copper, data=pbc2, ties="breslow")
royston(cfit)
}
\keyword{ survival }
