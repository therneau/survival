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
  
  An adjustment for D is based on the ratio
  r= (number of events)/(number of coefficients). For models which have
  sufficient sample size (r>20) the adjustment will be small.

  The Nagelkirke value is the Cox-Snell R-squared divided by a scaling
  constant. The two separate values are present in the result of
  \code{summary.coxph} as a 2 element vector \code{rsq}, and were listed as
  "Rsquare" and "max possible" in older versions of the print routine.
  (Since superseded in the default printout by the concordance.)
  The Nagelkirke estimate is not returned when \code{newdata} is present.
}
\value{a vector containing the value of D, the estimated standard error
  of D, and three or four pseudo R-squared values.
}
\references{
M. Goen and G. Heller, Concordance probability and discriminatory power
in proportional hazards regression.  Biometrika 92:965-970, 2005.

N. Nagelkerke, J. Oosting, J. and A. Hart, A simple test for goodness of fit
of Cox's proportional hazards model.  Biometrics 40:483-486, 1984.

P. Royston and W. Sauerbrei, A new measure of prognostic separation in
survival data.  Statistics in Medicine 23:723-748, 2004.
}
\examples{
# An example used in Royston and Sauerbrei
pbc2 <- na.omit(pbc)  # no missing values
cfit <- coxph(Surv(time, status==2) ~ age + log(bili) + edema + albumin +
                   stage + copper, data=pbc2, ties="breslow")
royston(cfit)
}
\keyword{ survival }
