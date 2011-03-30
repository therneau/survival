\name{pspline}
\alias{pspline}
\alias{psplineinverse}
\title{Smoothing splines using a pspline basis}
\usage{
pspline(x, df=4, theta, nterm=2.5 * df, degree=3, eps=0.1, method,
   Boundary.knots=range(x), intercept=FALSE, penalty=TRUE,  ...)

psplineinverse(x)}
\arguments{
  \item{x}{for psline: a covariate vector.
    The function does not apply to factor variables.
    For psplineinverse x will be the result of a pspline call.}
  \item{df}{the desired degrees of freedom.
    One of the arguments \code{df} or \code{theta}' must be given, but not both.
    If \code{df=0}, then the AIC = (loglik -df) is used to choose an
    "optimal" degrees of freedom.  If AIC is chosen, then an optional
    argument `caic=T' can be used to specify the corrected AIC of
    Hurvich et. al.
  }
  \item{theta}{roughness penalty for the fit.
    It is a monotone function of the degrees of freedom, with theta=1
    corresponding to a linear fit and theta=0 to an unconstrained fit
    of nterm degrees of freedom.
  }
  \item{nterm}{ number of splines in the basis }
  \item{degree}{ degree of splines }
  \item{eps}{accuracy for \code{df} }
  \item{method}{the method for choosing the tuning parameter \code{theta}.
    If theta is given, then 'fixed' is assumed.
    If the degrees of freedom is given, then 'df' is assumed.
    If method='aic' then the degrees of freedom is chosen automatically using
    Akaike's information criterion.}
  \item{\dots}{optional arguments to the control function}
  \item{Boundary.knots}{the spline is linear beyond the boundary knots.
    These default to the range of the data.}
  \item{intercept}{if TRUE, the basis functions include the intercept.}
  \item{penalty}{if FALSE a large number of attributes having to do with
    penalized fits are excluded.  Most useful for exploring the code so
    as to return a matrix with few added attributes.}
  }
\description{
    Specifies a penalised spline basis for the predictor.
    This is done by fitting a comparatively small set of splines and
    penalising the integrated second derivative.
    Traditional smoothing splines use one basis per observation, but several
    authors have pointed out that the final results of the fit are 
    indistinguishable for any number of basis functions greater than about 
    2-3 times the degrees of freedom. 
    Eilers and Marx point out that if the basis functions are evenly spaced,
    this leads to significant computational simplifications.
}

\value{
  Object of class \code{pspline, coxph.penalty} containing the spline basis,
 with the appropriate attributes to be
 recognized as a penalized term by the coxph or survreg functions.

 For psplineinverse the original x vector is reconstructed.
}

\seealso{\code{\link{coxph}},\code{\link{survreg}},\code{\link{ridge}},
  \code{\link{frailty}}
}
\references{
Eilers, Paul H. and Marx, Brian D. (1996).
Flexible smoothing with B-splines and penalties.
    Statistical Science, 11, 89-121.

Hurvich, C.M. and Simonoff, J.S. and Tsai, Chih-Ling (1998).
Smoothing parameter selection in nonparametric regression using
        an improved Akaike information criterion,
JRSSB, volume 60, 271--293.
}    
\examples{
lfit6 <- survreg(Surv(time, status)~pspline(age, df=2), cancer)
plot(cancer$age, predict(lfit6), xlab='Age', ylab="Spline prediction")
title("Cancer Data")
fit0 <- coxph(Surv(time, status) ~ ph.ecog + age, cancer)
fit1 <- coxph(Surv(time, status) ~ ph.ecog + pspline(age,3), cancer)
fit3 <- coxph(Surv(time, status) ~ ph.ecog + pspline(age,8), cancer)
fit0
fit1
fit3
}
\keyword{ survival}

