\name{survival-deprecated}
\alias{survival-deprecated}
\alias{survConcordance}
\alias{survConcordance.fit}
\title{Deprecated functions in package \pkg{survival}}
\description{
  These functions are temporarily retained for compatability with older programs,
  and may transition to defunct status.
}
\usage{
survConcordance(formula, data, weights, subset, na.action) # use concordance
survConcordance.fit(y, x, strata, weight)    # use concordancefit
}
\arguments{
  \item{formula}{
    a formula object, with the response on the left of a \code{~} operator, and 
    the terms on the right.  The response must be a survival object as 
    returned by the \code{Surv} function. 
  }
\item{data}{a data frame
}
\item{weights,subset,na.action}{as for \code{coxph}}
\item{x, y, strata, weight}{predictor, response, strata, and weight
  vectors for the direct call}
}
\seealso{
  \code{\link{Deprecated}}
}
\keyword{survival}