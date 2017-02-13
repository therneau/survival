\name{predict.survreg}
\alias{predict.survreg}
\alias{predict.survreg.penal}
\title{
  Predicted Values for a `survreg' Object
}
\description{
Predicted values for a \code{survreg} object 
}
\usage{
\method{predict}{survreg}(object, newdata,  
type=c("response", "link", "lp", "linear", "terms", "quantile",  
	"uquantile"),
 se.fit=FALSE, terms=NULL, p=c(0.1, 0.9), na.action=na.pass, ...)
}
\arguments{
\item{object}{
result of a model fit using the \code{survreg} function. 
}
\item{newdata}{
data for prediction.  If absent predictions are for the 
subjects used in the original fit. 
}
\item{type}{
the type of predicted value.  
This can be on the original scale of the data (response), 
the linear predictor (\code{"linear"}, with \code{"lp"} as an allowed abbreviation), 
a predicted quantile on the original scale of the data (\code{"quantile"}), 
a quantile on the linear predictor scale (\code{"uquantile"}),
or the matrix of terms for the linear predictor (\code{"terms"}).
At this time \code{"link"} and linear predictor (\code{"lp"}) are identical. 
}
\item{se.fit}{
if \code{TRUE}, include the standard errors of the prediction in the result. 
}
\item{terms}{
subset of terms.  The default for residual type \code{"terms"} is a matrix with 
one column for every term (excluding the intercept) in the model. 
}
\item{p}{
vector of percentiles.  This is used only for quantile predictions. 
}
\item{na.action}{
  applies only when the \code{newdata} argument is present, and defines
  the missing value action for the new data.  The default is to include
  all observations.}

\item{\dots}{for future methods}
}
\value{
a vector or matrix of predicted values. 
}
\references{
Escobar and Meeker (1992). Assessing influence in regression analysis with 
censored data. \emph{Biometrics,} 48, 507-528. 
}
\seealso{
\code{\link{survreg}},  \code{\link{residuals.survreg}}   
}
\examples{
# Draw figure 1 from Escobar and Meeker, 1992.
fit <- survreg(Surv(time,status) ~ age + I(age^2), data=stanford2, 
	dist='lognormal') 
with(stanford2, plot(age, time, xlab='Age', ylab='Days', 
	xlim=c(0,65), ylim=c(.1, 10^5), log='y', type='n'))
with(stanford2, points(age, time, pch=c(2,4)[status+1], cex=.7))
pred <- predict(fit, newdata=list(age=1:65), type='quantile', 
	         p=c(.1, .5, .9)) 
matlines(1:65, pred, lty=c(2,1,2), col=1) 

# Predicted Weibull survival curve for a lung cancer subject with
#  ECOG score of 2
lfit <- survreg(Surv(time, status) ~ ph.ecog, data=lung)
pct <- 1:98/100   # The 100th percentile of predicted survival is at +infinity
ptime <- predict(lfit, newdata=data.frame(ph.ecog=2), type='quantile',
                 p=pct, se=TRUE)
matplot(cbind(ptime$fit, ptime$fit + 2*ptime$se.fit,
                         ptime$fit - 2*ptime$se.fit)/30.5, 1-pct,
        xlab="Months", ylab="Survival", type='l', lty=c(1,2,2), col=1)
}
\keyword{survival}
