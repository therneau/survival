\name{residuals.survreg}
\alias{residuals.survreg}
\alias{residuals.survreg.penal}
\title{Compute Residuals for `survreg' Objects}
\description{
This is a method for the function \code{\link{residuals}} for objects
inheriting from class \code{survreg}.  
}
\usage{
\method{residuals}{survreg}(object, type=c("response", "deviance","dfbeta","dfbetas",
"working","ldcase","ldresp","ldshape", "matrix"), rsigma=TRUE,
collapse=FALSE, weighted=FALSE, ...)
}
\arguments{
\item{object}{
an object inheriting from class \code{survreg}.
}
\item{type}{
type of residuals, with choices of \code{"response"}, \code{"deviance"},
\code{"dfbeta"}, \code{"dfbetas"}, \code{"working"}, \code{"ldcase"}, \code{"lsresp"},
\code{"ldshape"}, and \code{"matrix"}.  See the LaTeX documentation
(\code{survival/doc/survival.ps.gz}) for more
detail.
}
\item{rsigma}{
include the scale parameters in the variance matrix, when doing computations.
(I can think of no good reason not to).
}
\item{collapse}{
optional vector of subject groups.  If given, this must be of the same
length as the residuals, and causes the result to be per group residuals.
}
\item{weighted}{
give weighted residuals?  Normally residuals are unweighted.
}\item{...}{other unused arguments}}
\value{
A vector or matrix of residuals is returned.
Response residuals are on the scale of the original data,
working residuals are on the scale of the linear predictor,
and deviance residuals are on log-likelihood scale.
The dfbeta residuals are a matrix, where the ith row gives the
approximate change in the coefficients due to the addition of subject i.
The dfbetas matrix contains the dfbeta residuals, with each column
scaled by the standard deviation of that coefficient.


The matrix type produces a matrix based on derivatives of the log-likelihood
function.  Let \eqn{L} be the log-likelihood, \eqn{p} be the linear predictor \eqn{X\beta}{X \%*\% coef},
and \eqn{s} be \eqn{\log(\sigma)}.  Then the 6 columns of the matrix are
\eqn{L}, \eqn{dL/dp},\eqn{\partial^2L/\partial p^2}{ddL/(dp dp)},
\eqn{dL/ds}, \eqn{\partial^2L/\partial s^2}{ddL/(ds ds)}  and
\eqn{\partial^2L/\partial p\partial s}{ddL/(dp ds)}.  Diagnostics based on these quantities
are discussed in an article by Escobar and Meeker.
The main ones are the likelihood displacement residuals for perturbation
of a case weight (\code{ldcase}), the response value (\code{ldresp}), and the \code{shape}.
}
\references{
Escobar, L. A. and Meeker, W. Q. (1992).
Assessing influence in regression analysis with censored data.
\emph{Biometrics}
\bold{48}, 507-528.
}
\seealso{\code{\link{predict.survreg}}}
\examples{
fit <- survreg(Surv(time,status) ~x, aml)
rr  <- residuals(fit, type='matrix')
}
\keyword{survival}
% Converted by Sd2Rd version 0.3-2.
