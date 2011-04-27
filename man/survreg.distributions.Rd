\name{survreg.distributions}
\alias{survreg.distributions}
\title{Parametric Survival Distributions}
\usage{
survreg.distributions
}
\description{ List of distributions for accelerated failure models. These are
  location-scale families for some transformation of time. The entry
  describes  the cdf \eqn{F} and density \eqn{f} of a canonical member of
  the family.
}

\format{
 
  There are two basic formats, the first defines a distribution de novo,
the second defines a new distribution in terms of an old one.
  \tabular{ll}{
    name:\tab name of distribution\cr
    variance:\tab function(parms) returning the variance (currently unused)\cr
    init(x,weights,...):\tab Function returning an initial\cr
    \tab estimate of the mean and variance \cr
    \tab (used for initial values in the iteration)\cr
    density(x,parms):\tab Function returning a matrix with columns \eqn{F},\cr
    \tab \eqn{1-F},\eqn{f},\eqn{f'/f},\eqn{f''/f}\cr
    quantile(p,parms):\tab Quantile function\cr
    scale:\tab Optional fixed value for the scale parameter\cr
    parms:\tab Vector of default values and names for any additional parameters\cr
    deviance(y,scale,parms):\tab Function returning the deviance for a\cr
    \tab saturated model; used only for deviance residuals.
  }
  and to define one distribution in terms of another
  \tabular{ll}{
    name:\tab name of distribution\cr
    dist:\tab name of parent distribution\cr
    trans:\tab transformation (eg log)\cr
    dtrans:\tab derivative of transformation\cr
    itrans:\tab inverse of transformation\cr
    scale:\tab Optional fixed value for scale parameter\cr
  }
}

\details{
There are four basic distributions:\code{extreme}, \code{gaussian},
\code{logistic} and \code{t}. The last three
are parametrised in the same way as the distributions already present in
\R. The extreme value cdf is
\deqn{F=1-e^{-e^t}.}

When the logarithm of survival time has one of the first three distributions
we obtain respectively \code{weibull}, \code{lognormal}, and
\code{loglogistic}. The location-scale parameterizaion of a Weibull
distribution found in \code{survreg} is not the same as the parameterization
of \code{\link{rweibull}}.

The other predefined distributions are defined in terms of these. The
\code{exponential} and \code{rayleigh} distributions are Weibull
distributions with fixed \code{scale} of 1 and 0.5 respectively, and
\code{loggaussian} is a synonym for \code{lognormal}.

For speed parts of the three most commonly used distributions
are hardcoded in C; for this reason the elements of \code{survreg.distributions}
with names of "Extreme value", "Logisitic" and "Gaussian" should not be
modified.  (The order of these in the list is not important, recognition
is by name.)
As an alternative to modifying \code{survreg.distributions}
a new distribution can be specified as a separate list.
This is the preferred method of addition and is illustrated below.
}
\seealso{\code{\link{survreg}}, \code{\link{pweibull}},
  \code{\link{pnorm}},\code{\link{plogis}}, \code{\link{pt}},
  \code{\link{survregDtest}}
}
\examples{
# time transformation
survreg(Surv(time, status) ~ ph.ecog + sex, dist='weibull', data=lung)
# change the transformation to work in years
# intercept changes by log(365), everything else stays the same
my.weibull <- survreg.distributions$weibull
my.weibull$trans <- function(y) log(y/365)
my.weibull$itrans <- function(y) 365*exp(y)
survreg(Surv(time, status) ~ ph.ecog + sex, lung, dist=my.weibull)

# Weibull parametrisation
y<-rweibull(1000, shape=2, scale=5)
survreg(Surv(y)~1, dist="weibull")
# survreg parameters are scale=1/shape, intercept=log(scale)

# Cauchy fit
mycauchy <- list(name='Cauchy',
                 init= function(x, weights, ...) 
                      c(median(x), mad(x)),
                 density= function(x, parms) {
                      temp <- 1/(1 + x^2)
                      cbind(.5 + atan(x)/pi, .5+ atan(-x)/pi,
                            temp/pi, -2 *x*temp, 2*temp*(4*x^2*temp -1))
                      },
                 quantile= function(p, parms) tan((p-.5)*pi),
                 deviance= function(...) stop('deviance residuals not defined')
                 )
survreg(Surv(log(time), status) ~ ph.ecog + sex, lung, dist=mycauchy)
}
\keyword{survival}
