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
\code{"ldshape"}, and \code{"matrix"}. 
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
\eqn{\partial^2L/\partial p\partial s}{ddL/(dp ds)}.
Diagnostics based on these quantities
are discussed in the book and article by Escobar and Meeker.
The main ones are the likelihood displacement residuals for perturbation
of a case weight (\code{ldcase}), the response value (\code{ldresp}),
and the \code{shape}.

For a transformed distribution such as the log-normal or Weibull,
matrix residuals are based on the log-likelihood of the transformed data
log(y). 
For a monotone function f the density of f(X) is the density of X
divided by the derivative of f (the Jacobian), so subtract log(derivative) from
each uncensored observation's loglik value in order to match the
\code{loglik} component of the result.  The other colums of the matrix
residual are unchanged by the transformation.
}
\references{
Escobar, L. A. and Meeker, W. Q. (1992).
Assessing influence in regression analysis with censored data.
\emph{Biometrics}
\bold{48}, 507-528.

Escobar, L. A. and Meeker, W. Q. (1998).
Statistical Methods for Reliablilty Data.  Wiley.
}

\seealso{\code{\link{predict.survreg}}}
\examples{
fit <- survreg(Surv(futime, death) ~ age + sex, mgus2)
summary(fit)   # age and sex are both important

rr  <- residuals(fit, type='matrix')
sum(rr[,1]) - with(mgus2, sum(log(futime[death==1]))) # loglik

plot(mgus2$age, rr[,2], col= (1+mgus2$death)) # ldresp
}
\keyword{survival}

