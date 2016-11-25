\name{cipoisson}
\alias{cipoisson}
\title{Confidence limits for the Poisson}
\usage{
cipoisson(k, time = 1, p = 0.95, method = c("likelihood", "anscombe"))
}
\arguments{
\item{k}{Number of successes}

\item{time}{Total time on trial}

\item{p}{Probability level for the (two-sided) interval}

\item{method}{The method for computing the interval.}
}
\value{a vector, matrix, or array.
  If both \code{k} and \code{time} are single values the result is a
  vector of length 2 containing the lower an upper limits.
  If either or both are vectors the result is a matrix with two columns.
  If \code{k} is a matrix or array, the result will be an array with one
  more dimension; in this case the dimensions and dimnames (if any) of
  \code{k} are preserved.
}
\details{
The likelihood method is based on equation 10.10 of Feller, which relates
poisson probabilities to tail area of the gamma distribution.
The Anscombe approximation is based on the fact that sqrt(k + 3/8) is
has a nearly constant variance of 1/4, along with a continuity correction.
}
\examples{
cipoisson(4) # 95\\\% confidence limit 
# lower    upper  
# 1.089865 10.24153 
ppois(4, 10.24153)     #chance of seeing 4 or fewer events with large rate  
# [1] 0.02500096 
1-ppois(3, 1.08986)    #chance of seeing 4 or more, with a small rate 
# [1] 0.02499961

}
\references{
W.F. Feller (1950). An Introduction to Probability Theory and its
Applications, Volume 1, Chapter 6, Wiley.

F.J. Anscombe (1949). Transformations of Poisson, binomial and
negative-binomial data. Biometrika, vol 35, p246-254.
}
\seealso{
\code{\link[stats]{ppois}}, \code{\link[stats]{qpois}}
}

