\name{nwtco}
\alias{nwtco}
\docType{data}
\title{Data from the National Wilm's Tumor Study}
\description{
  Measurement error example. Tumor histology predicts
  survival, but prediction is stronger with central lab histology than
  with the local institution determination.
}
\usage{nwtco
}
\format{
  A data frame with 4028 observations on the following 9 variables.
  \describe{
    \item{\code{seqno}}{id number}
    \item{\code{instit}}{Histology from local institution}
    \item{\code{histol}}{Histology from central lab}
    \item{\code{stage}}{Disease stage}
    \item{\code{study}}{study}
    \item{\code{rel}}{indicator for relapse}
    \item{\code{edrel}}{time to relapse}
    \item{\code{age}}{age in months}
    \item{\code{in.subcohort}}{Included in the subcohort for the example in the
    paper}
  }
}
\references{
  NE Breslow and N Chatterjee (1999),
  Design and analysis of two-phase studies with binary outcome applied
  to Wilms tumour prognosis.
  \emph{Applied Statistics} \bold{48}, 457--68.
}
\note{
  The usage above assumes that the survival package is attached via
  \code{library(survival)} or via another package that depends on
  survival, e.g., mstate.  If not, then either \code{survival::nwtco}
  can be used or \code{data(nwtco, package="survival")}.
}
\examples{
with(nwtco, table(instit,histol))
anova(coxph(Surv(edrel,rel)~histol+instit,data=nwtco))
anova(coxph(Surv(edrel,rel)~instit+histol,data=nwtco))
}
\keyword{datasets}
