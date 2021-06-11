\name{royston}
\alias{royston}
\title{Compute Royston's D for a Cox model}
\description{
  Compute the D statistic proposed by Royston and Sauerbrei along with
  several pseudo- R square values.
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
  These values are called pseudo R-squared since they involve only
  the linear predictor, and not the outcome.
  \code{R.D} is the value that corresponsds the Royston and Sauerbrei
  \eqn{D}{D} statistic.  \code{R.KO} is the value proposed by Kent and
  O'Quigley, \code{R.N} is the value proposed by Nagelkerke, and
  \code{C.GH} corresponds to Goen and Heller's concordance measure.
  See the vignette for further details.
  
  An adjustment for D is based on the ratio
  r= (number of events)/(number of coefficients). For models which have
  sufficient sample size (r>20) the adjustment will be small. 
}
\value{a vector containing the value of D, the estimated standard error
  of D, and four pseudo R-squared values.
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
