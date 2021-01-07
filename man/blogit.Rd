\name{blogit}
\alias{blogit}
\alias{bcloglog}
\alias{bprobit}
\title{
Bounded link functions
}
\description{
Alternate link functions that impose bounds on the input of their link function
}
\usage{
blogit(edge = 0.05)
bprobit(edge= 0.05)
bcloglog(edge=.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{edge}{input is mapped to the range of (edge, 1-edge)}
}
\details{
  When using survival psuedovalues for regression, the raw data can be
  outside the range (0,1), yet we want to restrict the predicted values
  to lie within that range.  A natural way to deal with this is to use
  \code{glm} with \code{family = gaussian(link= "logit")}.
  Butthis will fail.
  The reason is that the \code{family} object has a component
  \code{linkfun} that does not accept values outside of (0,1).
  
  This function is only used to create initial values for the iteration
  step, however. Mapping the offending input argument into the range
  of (egde, 1-edge) before computing the link results in starting
  estimates that are good enough.  The final result of the fit will be
  no different than if explicit starting estimates were given using the
  \code{etastart} or \code{mustart} arguments.
  These functions create copies of the logit, probit, and complimentary
  log-log families that differ from the standard ones only in this
  use of a bounded input argument, and are called a "bounded logit" =
  \code{blogit}, etc.
  
  A \code{family} object used by the \code{glm} function has a component
  \code{linkfun}, which is used only to establish initial values for the
  iteration.
  When using survival psuedovalues for regression, the raw data can be
  outside the range (0,1), yet we want to restrict the predicted values
  to lie within that range, which suggests 
}
\value{a \code{family} object of the same form as \code{make.family}.
}
\seealso{\code{\link{stats}{make.family}}}
\examples{
py <- pseudo(survfit(Surv(time, status) ~1, lung), time=730) #2 year survival
range(py)
pfit <- glm(py ~ ph.ecog, data=lung, family=gaussian(link=blogit()))
# For each +1 change in performance score, the odds of 2 year survival
#  are multiplied by 1/2  (exp of the coefficient) 
}
\keyword{survival}