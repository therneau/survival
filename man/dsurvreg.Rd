\name{dsurvreg}
\alias{dsurvreg}
\alias{psurvreg}
\alias{qsurvreg}
\title{
  Distributions available in survreg.
}
\description{
Density, cumulative probability, and quantiles for the set of distributions
supported by the \code{survreg} function.
}
\usage{
dsurvreg(x, mean, scale=1, distribution='weibull', parms)
psurvreg(q, mean, scale=1, distribution='weibull', parms)
qsurvreg(p, mean, scale=1, distribution='weibull', parms)
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
\item{mean}{
vector of means. 
This is replicated to be the same length as \code{p} 
or \code{q}.
}
\item{scale}{
vector of (positive) scale factors.
This is replicated to be the same length as \code{p} 
or \code{q}.
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

The \code{mean} and \code{scale}
values are as they would be for \code{survreg}.
In particular, if
the distribution is one that involves a transformation, then they are the
mean and scale of the transformed distribution.  
For example, the Weibull distribution is fit using the
Extreme value distribution along with a log transformation.
Letting \eqn{F(t) = 1 - \exp[-(at)^p]}{F(t) = 1 - exp(-(at)^p)}
be the cumulative distribution of the
Weibull, the mean corresponds to \eqn{-\log(a)}{-log(a)} and the scale
to \eqn{1/p} 
(Kalbfleish and Prentice, section 2.2.2).
}
\section{References}{
Kalbfleish, J. D. and Prentice, R. L. (1970). 
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
}
\keyword{distribution}
% docclass is function
% Converted by Sd2Rd version 37351.
