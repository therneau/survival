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
se.fit=FALSE, na.action=na.pass, terms=names(object$assign), collapse, ...)
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
  all observations.}
\item{terms}{
if type="terms", this argument can be used to specify which terms should be 
included; the default is all. 
}
\item{collapse}{
optional vector of subject identifiers.  
If specified, the output will contain one entry per subject rather than one 
entry per observation. 
}
\item{\dots}{For future methods}
}
\value{
a vector or matrix of predictions, or a list containing the predictions 
(element "fit") and their standard errors (element "se.fit") if the se.fit 
option is TRUE. 
}

\seealso{
\code{\link{predict}},\code{\link{coxph}},\code{\link{termplot}} }

\examples{
fit <- coxph(Surv(time, status) ~ age + ph.ecog + strata(inst), lung) 
mresid <- lung$status - predict(fit, type='expected') #Martingale resid 
predict(fit,type="lp")
predict(fit,type="risk")
predict(fit,type="expected")
predict(fit,type="terms")
predict(fit,type="lp",se.fit=TRUE)
predict(fit,type="risk",se.fit=TRUE)
predict(fit,type="expected",se.fit=TRUE)
predict(fit,type="terms",se.fit=TRUE)
}
\keyword{survival}
