\name{tobin}
\alias{tobin}
\docType{data}
\title{Tobin's Tobit data}
\description{
  Economists fit a parametric censored data model called the
  \sQuote{tobit}.  These data are from Tobin's original paper. 
}
\usage{tobin}
\format{
  A data frame with 20 observations on the following 3 variables.
  \describe{
    \item{durable}{Durable goods purchase}
    \item{age}{Age in years}
    \item{quant}{Liquidity ratio (x 1000)}
  }
}
\source{
  J Tobin (1958),
  Estimation of relationships for limited dependent variables.
  \emph{Econometrica} \bold{26}, 24--36.
}
\examples{
tfit <- survreg(Surv(durable, durable>0, type='left') ~age + quant,
                data=tobin, dist='gaussian')

predict(tfit,type="response")

}
\keyword{datasets}
