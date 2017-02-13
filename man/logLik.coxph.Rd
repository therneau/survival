\name{logLik.coxph}
\alias{logLik.coxph}
\alias{logLik.survreg}
\title{logLik method for a Cox model}
\description{The logLik function for survival models}
\usage{
\method{logLik}{coxph}(object, ...)
\method{logLik}{survreg}(object, ...)
}
\arguments{
  \item{object}{the result of a \code{coxph} or \code{survreg} fit}
  \item{\dots}{optional arguments for other instances of the method}
}
\details{
 The logLik function is used by summary functions in R such as
 \code{AIC}.
 For a Cox model, this method returns the partial likelihood.
 The number of degrees of freedom (df) used by the fit and the effective
 number of observations (nobs) are added as attributes.
 Per Raftery and others, the effective number of observations is the
 taken to be the number of events in the data set.

 For a \code{survreg} model the proper value for the effective number
 of observations is still an open question (at least to this author).
 For right censored data the approach of \code{logLik.coxph} is the
 possible the most sensible, but for interval censored observations
 the result is unclear.  The code currently does not add a \emph{nobs}
 attribute.
 }
 \value{an object of class \code{logLik}}
 
 \references{
   Robert E. Kass and Adrian E. Raftery (1995). "Bayes Factors". J.
   American Statistical Assoc. 90 (430): 791.

   Raftery A.E. (1995), "Bayesian Model Selection in Social Research",
   Sociological methodology, 111-196.
 }
 \seealso{\code{\link{logLik}}}
 \author{Terry Therneau}

\keyword{ survival}
