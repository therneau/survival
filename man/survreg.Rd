\name{survreg}
\alias{survreg}
\alias{model.frame.survreg}
\alias{labels.survreg}
\alias{print.survreg.penal}
\alias{print.summary.survreg}
\alias{survReg}
\alias{anova.survreg}
\alias{vcov.survreg}
\alias{anova.survreglist}
\title{
  Regression for a Parametric Survival Model 
}
\description{
Fit a parametric survival regression model.
These are location-scale models for an arbitrary transform of the time
variable; the most common cases use a log transformation, leading to
accelerated failure time models.
}
\usage{
survreg(formula, data, weights, subset, 
        na.action, dist="weibull", init=NULL, scale=0, 
        control,parms=NULL,model=FALSE, x=FALSE,
        y=TRUE, robust=FALSE, score=FALSE, \dots)

}
\arguments{
\item{formula}{
a formula expression as for other regression models. 
The response is usually a survival object as returned by the \code{Surv} function. 
See the documentation for \code{Surv}, \code{lm} and \code{formula} for details. 
}
\item{data}{
a data frame in which to interpret the variables named in 
the \code{formula}, \code{weights} or the \code{subset} arguments. 
}
\item{weights}{optional vector of case weights}
\item{subset}{
  subset of the observations to be used in the fit
}
\item{na.action}{
a missing-data filter function, applied to the model.frame, after any 
\code{subset} argument has been used.  Default is \code{options()\$na.action}. 
}
\item{dist}{
assumed distribution for y variable. 
If the argument is a character string, then it is assumed to name an
element from \code{\link{survreg.distributions}}. These include
\code{"weibull"}, \code{"exponential"}, \code{"gaussian"},
\code{"logistic"},\code{"lognormal"} and \code{"loglogistic"}.
Otherwise, it is assumed to be a user defined list conforming to the
format described in \code{\link{survreg.distributions}}.
}
\item{parms}{
a list of fixed parameters.  For the t-distribution for instance this is
the degrees of freedom; most of the distributions have no parameters.
}
\item{init}{
optional vector of initial values for the parameters.
}
\item{scale}{
optional fixed value for the scale.  If set to <=0 then the scale is
estimated.
}
\item{control}{
a list of control values, in the format produced by
\code{\link{survreg.control}}. The default value is \code{survreg.control()}
}
\item{model,x,y}{
flags to control what is returned.  If any of these is true, then the
model frame, the model matrix, and/or the vector of response times will be
returned as components of the final result, with the same names as the
flag arguments.}
\item{score}{return the score vector. (This is expected to be zero upon
successful convergence.)
}
\item{robust}{Use robust 'sandwich' standard errors, based on
  independence of individuals if there is no \code{cluster()} term in
  the formula, based on independence of clusters if there is.}
\item{\dots}{
other arguments which will be passed to \code{survreg.control}.
}}
\value{
an object of class \code{survreg} is returned.
}
\seealso{
\code{\link{survreg.object}}, \code{\link{survreg.distributions}},
\code{\link{pspline}}, \code{\link{frailty}}, \code{\link{ridge}}
}
\examples{
# Fit an exponential model: the two fits are the same
survreg(Surv(futime, fustat) ~ ecog.ps + rx, ovarian, dist='weibull',
                                    scale=1)
survreg(Surv(futime, fustat) ~ ecog.ps + rx, ovarian,
        dist="exponential")

#
# A model with different baseline survival shapes for two groups, i.e.,
#   two different scales
survreg(Surv(time, status) ~ ph.ecog + age + strata(sex), lung)

# There are multiple ways to parameterize a Weibull distribution. The survreg 
# function imbeds it in a general location-scale familiy, which is a 
# different parameterization than the rweibull function, and often leads
# to confusion.
#   survreg's scale  =    1/(rweibull shape)
#   survreg's intercept = log(rweibull scale)
#   For the log-likelihood all parameterizations lead to the same value.
y <- rweibull(1000, shape=2, scale=5)
survreg(Surv(y)~1, dist="weibull")

# Economists fit a model called `tobit regression', which is a standard
# linear regression with Gaussian errors, and left censored data.
tobinfit <- survreg(Surv(durable, durable>0, type='left') ~ age + quant,
	            data=tobin, dist='gaussian')
}
\keyword{survival}
