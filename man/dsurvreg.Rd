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
\item{mean}{vector of linear predictors for the model.  
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
since almost none of these distributions are symmetric it will not
actually be a mean, but corresponds instead to the linear predictor of
a fitted model.
Translation to the usual parameterization found in a textbook is not
always obvious.
For example, the Weibull distribution is fit using the
Extreme value distribution along with a log transformation.
Letting \eqn{F(t) = 1 - \exp[-(at)^p]}{F(t) = 1 - exp(-(at)^p)}
be the cumulative distribution of the
Weibull using a standard parameterization in terms of
\eqn{a} and \eqn{p},
the survreg location corresponds to \eqn{-\log(a)}{-log(a)} and the scale
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
