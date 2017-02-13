\name{survreg.control}
\alias{survreg.control}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Package options for survreg and coxph}
\description{
This functions checks and packages the fitting options for 
\code{\link{survreg}}
}
\usage{
survreg.control(maxiter=30, rel.tolerance=1e-09, 
toler.chol=1e-10, iter.max, debug=0, outer.max=10)

}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{maxiter}{maximum number of iterations }
  \item{rel.tolerance}{relative tolerance to declare convergence }
  \item{toler.chol}{Tolerance to declare Cholesky decomposition singular}
  \item{iter.max}{same as \code{maxiter}}
  \item{debug}{print debugging information}
  \item{outer.max}{maximum number of outer iterations for choosing
    penalty parameters}
}
\value{
 A list with the same elements as the input
}

\seealso{ \code{\link{survreg}}}
\keyword{survival}

