\name{coxph.wtest}
\alias{coxph.wtest}
\title{Compute a quadratic form}
\description{
  This function is used internally by several survival routines.  It
  computes a simple quadratic form, while properly dealing with missings.
}
\usage{
coxph.wtest(var, b, toler.chol = 1e-09)
}
\arguments{
  \item{var}{variance matrix}
  \item{b}{vector}
  \item{toler.chol}{tolerance for the internal choelsky decomposition}
}
\details{
  Compute b' V-inverse b.  Equivalent to sum(b * solve(V,b)), except for
  the case of redundant covariates in the original model, which lead to
  NA values in V and b.
}
\value{a real number}
\author{Terry Therneau}
\keyword{ survival }

