\name{ridge}
\alias{ridge}
\title{ Ridge regression}
\usage{
ridge(..., theta, df=nvar/2, eps=0.1, scale=TRUE)
}
\arguments{
 \item{\dots}{predictors to be ridged }
 \item{theta}{penalty is \code{theta}/2 time sum of squared coefficients }
 \item{df}{Approximate degrees of freedom }
 \item{eps}{ Accuracy required for \code{df} }
 \item{scale}{ Scale variables before applying penalty? }
}
\description{
 When used in a \link{coxph} or \link{survreg} model formula,
specifies a ridge regression term.  The likelihood is penalised by
\code{theta}/2 time the sum of squared coefficients. If \code{scale=T}
the penalty is calculated for coefficients based on rescaling the
predictors to have unit variance. If \code{df} is specified then \code{theta} is chosen based on an approximate degrees of freedom.
 } 

\value{
 An object of class \code{coxph.penalty} containing the data and
control functions.  } 

\references{
Gray (1992) "Flexible methods of analysing survival data using splines, with applications to breast cancer prognosis" JASA 87:942--951
 }

\seealso{ \code{\link{coxph}},\code{\link{survreg}},\code{\link{pspline}},\code{\link{frailty}} }

\examples{

coxph(Surv(futime, fustat) ~ rx + ridge(age, ecog.ps, theta=1),
	      ovarian)

lfit0 <- survreg(Surv(time, status) ~1, cancer)
lfit1 <- survreg(Surv(time, status) ~ age + ridge(ph.ecog, theta=5), cancer)
lfit2 <- survreg(Surv(time, status) ~ sex + ridge(age, ph.ecog, theta=1), cancer)
lfit3 <- survreg(Surv(time, status) ~ sex + age + ph.ecog, cancer)

}
\keyword{survival }%-- one or more ...
