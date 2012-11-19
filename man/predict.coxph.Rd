\name{predict.coxph}
\alias{predict.coxph}
\alias{predict.coxph.penal}
\title{
  Predictions for a Cox model 
}
\description{
  Compute fitted values and regression terms for a model fitted by
  \code{\link{coxph}}
}
\usage{
\method{predict}{coxph}(object, newdata,
type=c("lp", "risk", "expected", "terms"),
se.fit=FALSE, na.action=na.pass, terms=names(object$assign), collapse,
reference=c("strata", "sample"),  ...)
}
\arguments{
\item{object}{
the results of a coxph fit. 
}
\item{newdata}{
  Optional new data at which to do predictions.
  If absent predictions are for the 
data frame used in the original fit. 
}
\item{type}{
the type of predicted value. 
Choices are the linear predictor (\code{"lp"}), the risk score exp(lp)
(\code{"risk"}), 
the expected number of events given the covariates and follow-up time 
(\code{"expected"}), and the terms of the linear predictor (\code{"terms"}). 
}
\item{se.fit}{
if TRUE, pointwise standard errors are produced for the predictions. 
}
\item{na.action}{
  applies only when the \code{newdata} argument is present, and defines
  the missing value action for the new data.  The default is to include
  all observations.
  When there is no newdata, then the behavior of missing is dictated by
  the na.action option of the original fit.}
\item{terms}{
if type="terms", this argument can be used to specify which terms should be 
included; the default is all. 
}
\item{collapse}{
optional vector of subject identifiers.  
If specified, the output will contain one entry per subject rather than one 
entry per observation. 
}
\item{reference}{reference for centering predictions, see details below}
\item{\dots}{For future methods}
}
\value{
a vector or matrix of predictions, or a list containing the predictions 
(element "fit") and their standard errors (element "se.fit") if the se.fit 
option is TRUE. 
}
\details{
  The Cox model is a \emph{relative} risk model; predictions
of type "linear predictor", "risk", and "terms" are all
relative to the sample from which they came.  By default, the reference
value for each of these is the mean covariate within strata.  The
primary underlying
reason is statistical: a Cox model only predicts relative risks
between pairs of subjects within the same strata, and hence the addition
of a constant to any covariate, either overall or only within a
particular stratum, has no effect on the fitted results.
Using the \code{reference="strata"} option causes this to be true for
predictions as well.

When the results of \code{predict} are used in further calculations it
may be desirable to use a fixed reference level.
Use of \code{reference="sample"} will use the overall means, and agrees
with the \code{linear.predictors} component of the coxph object (which
uses the overall mean for backwards compatability with older code).
Predictions of \code{type="terms"} are almost invariably passed
forward to further calculation, so for these we default to using
the sample as the reference.

Predictions of type "expected" incorporate the baseline hazard and are
thus absolute instead of relative; the
\code{reference} option has no effect on these.

Models that contain a \code{frailty} term are a special case: due
to the technical difficulty, when there is a \code{newdata} argument the
predictions will always be for a random effect of zero.  
}

\seealso{
\code{\link{predict}},\code{\link{coxph}},\code{\link{termplot}} }

\examples{
options(na.action=na.exclude) # retain NA in predictions
fit <- coxph(Surv(time, status) ~ age + ph.ecog + strata(inst), lung)
#lung data set has status coded as 1/2
mresid <- (lung$status-1) - predict(fit, type='expected') #Martingale resid 
predict(fit,type="lp")
predict(fit,type="expected")
predict(fit,type="risk",se.fit=TRUE)
predict(fit,type="terms",se.fit=TRUE)
}
\keyword{survival}
