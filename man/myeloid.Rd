\name{myeloid}
\alias{myeloid}
\docType{data}
\title{Acute myeloid leukemia}
\description{
  This simulated data set is based on a trial in acute myeloid
  leukemia.
}
\format{
  A data frame with 646 observations on the following 9 variables.
  \describe{
    \item{\code{id}}{subject identifier, 1-646}
    \item{\code{trt}}{treatment arm A or B}
    \item{\code{sex}}{f=female, m=male}
    \item{\code{futime}}{time to death or last follow-up}
    \item{\code{death}}{1 if \code{futime} is a death, 0 for censoring}
    \item{\code{txtime}}{time to hematropetic stem cell transplant}
    \item{\code{crtime}}{time to complete response}
    \item{\code{rltime}}{time to relapse of disease}
  }
}
\details{
  This data set is used to illustrate multi-state survival curves.
  The correlation between within-subject event times strongly resembles
  that from an actual trial, but none of the actual data values are
  from that source.
}
\examples{
coxph(Surv(futime, death) ~ trt, data=myeloid)
# See the mstate vignette for a more complete analysis
}
\keyword{datasets}
\references{
  Le-Rademacher JG, Peterson RA, Therneau TM, Sanford BL, Stone RM,
  Mandrekar SJ.  Application of multi-state models in cancer clinical trials. 
  Clin Trials. 2018 Oct; 15 (5):489-498
  }