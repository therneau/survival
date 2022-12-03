\name{dsurvreg}
\alias{dsurvreg}
\alias{psurvreg}
\alias{qsurvreg}
\alias{rsurvreg}
\title{
  Distributions available in survreg.
}
\description{
  Density, cumulative distribution function, quantile function and random
  generation for the set of distributions
  supported by the \code{survreg} function.
}
\usage{
dsurvreg(x, mean, scale=1, distribution='weibull', parms)
psurvreg(q, mean, scale=1, distribution='weibull', parms)
qsurvreg(p, mean, scale=1, distribution='weibull', parms)
rsurvreg(n, mean, scale=1, distribution='weibull', parms)
}
\arguments{
\item{x}{
vector of quantiles. 
Missing values (\code{NA}s) are allowed. 
}
\item{q}{
vector of quantiles. 
Missing values (\code{NA}s) are allowed. 
}
\item{p}{
vector of probabilities. 
Missing values (\code{NA}s) are allowed. 
}
\item{n}{number of random deviates to produce}
\item{mean}{vector of location (linear predictor) parameters for the model.  
This is replicated to be the same length as \code{p}, \code{q}
or \code{n}.  
}
\item{scale}{
vector of (positive) scale factors.
This is replicated to be the same length as \code{p},
\code{q} or \code{n}.
}
\item{distribution}{
character string giving the name of the distribution.  This must be one
of the elements of \code{survreg.distributions}
}
\item{parms}{
optional parameters, if any, of the distribution.  For the t-distribution
this is the degrees of freedom.
}
}
\value{
density (\code{dsurvreg}), 
probability (\code{psurvreg}), 
quantile (\code{qsurvreg}), or 
for the requested distribution with mean and scale
parameters \code{mean} and 
\code{sd}. 
}
\details{
Elements of \code{q} or 
\code{p} that are missing will cause the corresponding 
elements of the result to be missing. 

The \code{location} and \code{scale}
values are as they would be for \code{survreg}.
The label "mean" was an unfortunate choice (made in mimicry of qnorm);
a more correct label would be "linear predictor".  
Since almost none of these distributions are symmetric the location
parameter is not actually a mean.

The \code{survreg} routines use the parameterization found in chapter
2 of Kalbfleisch and Prentice. 
Translation to the usual parameterization found in a textbook is not
always obvious.
For example, the Weibull distribution has cumulative distribution
function
\eqn{F(t) = 1 - e^{-(\lambda t)^p}}{F(t)}{1 - exp((-lambda t)^p)}.
The actual fit uses the fact that \eqn{\log(t)}{log(t)} has an extreme
value distribution, with location and scale of
\eqn{\alpha, \sigma}{a,s}, which are the location and scale parameters
reported by the \code{survreg} function.
The parameters are related by \eqn{\sigma= 1/p}{s=1/p} and
\eqn{\alpha = -\log(\lambda}{a= =log(lambda)}.
The \code{stats::dweibull} routine is parameterized in terms of
shape and scale parameters which correspond to \eqn{p}{p} and
\eqn{1/\lambda}{1/lambda} in the K and P notation.
Combining these we see that shape = \eqn{1/\sigma}{1/s} and
scale = \eqn{\exp{alpha}}{exp(a)}.
}
\section{References}{
Kalbfleisch, J. D. and Prentice, R. L. (1970).
\emph{The Statistical Analysis of Failure Time Data}
Wiley, New York. 
}
\seealso{
\code{\link{survreg}}, 
\code{\link{Normal}}
}
\examples{
# List of distributions available
names(survreg.distributions)
\dontrun{
 [1] "extreme"     "logistic"    "gaussian"    "weibull"     "exponential"
 [6] "rayleigh"    "loggaussian" "lognormal"   "loglogistic" "t"          
}
# Compare results
all.equal(dsurvreg(1:10, 2, 5, dist='lognormal'), dlnorm(1:10, 2, 5))

# Hazard function for a Weibull distribution
x   <- seq(.1, 3, length=30)
haz <- dsurvreg(x, 2, 3)/ (1-psurvreg(x, 2, 3))
\dontrun{
plot(x, haz, log='xy', ylab="Hazard") #line with slope (1/scale -1)
}

# Estimated CDF of a simple Weibull
fit <- survreg(Surv(time, status) ~ 1, data=lung)
pp <- 1:99/100  
q1 <- qsurvreg(pp, coef(fit), fit$scale)
q2 <- qweibull(pp, shape= 1/fit$scale, scale= exp(coef(fit)))
all.equal(q1, q2)
\dontrun{
plot(q1, pp, type='l', xlab="Months", ylab="CDF")
}
# per the help page for dweibull, the mean is scale * gamma(1 + 1/shape)
c(mean = exp(coef(fit))* gamma(1 + fit$scale))

}
\references{
  Kalbfleisch, J. D. and Prentice, R. L., The statistical analysis of
  failure time data, Wiley, 2002.
}
\keyword{distribution}
