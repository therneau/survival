\name{vcov.coxph}
\alias{vcov.coxph}
\alias{vcov.survreg}
\title{Variance-covariance matrix}
\description{Extract and return the variance-covariance matrix.}
\usage{
\method{vcov}{coxph}(object, complete=TRUE, ...)
\method{vcov}{survreg}(object, complete=TRUE, ...)
}
\arguments{
  \item{object}{a fitted model object}
  \item{complete}{logical indicating if the full variance-covariance
    matrix should be returned.  This has an effect only for an
    over-determined fit where some of the coefficients are undefined,
    and \code{coef(object)} contains corresponding NA values.
    If \code{complete=TRUE} the returned matrix will have row/column for
    each coefficient, if FALSE it will contain rows/columns
    corresponding to the non-missing coefficients.
    The coef() function has a simpilar \code{complete} argument.
  }
  \item{\ldots}{additional arguments for method functions}
}
\value{a matrix}
\details{
  For the \code{coxph} and \code{survreg} functions the returned matrix
is a particular generalized inverse:  the row and column corresponding
to any NA coefficients will be zero.  This is a side effect of the
generalized cholesky decomposion used in the unerlying compuatation.
}
\keyword{survival}

