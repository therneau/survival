\name{myeloid}
\alias{myeloid}
\docType{data}
\title{Acute myeloid leukemia}
\description{
  This simulated data set is based on a trial in acute myeloid
  leukemia.
}
\usage{myeloid
data(cancer, package="survival")
}
\format{
  A data frame with 646 observations on the following 9 variables.
  \describe{
    \item{\code{id}}{subject identifier, 1-646}
    \item{\code{trt}}{treatment arm A or B}
    \item{\code{sex}}{f=female, m=male}
    \item{\code{flt3}}{mutations of the FLT3 gene, a factor with levels
      of A, B, C}
    \item{\code{futime}}{time to death or last follow-up}
    \item{\code{death}}{1 if \code{futime} is a death, 0 for censoring}
    \item{\code{txtime}}{time to hematropetic stem cell transplant}
    \item{\code{crtime}}{time to complete response}
    \item{\code{rltime}}{time to relapse of disease}
  }
}
\details{
  This data set is used to illustrate multi-state survival curves.
  It is based on the actual study in the reference below.
  A subset of subjects was de-identifed, reordered, and then all of the
  time values randomly perturbed.

  Mutations in the FLT3 domain occur in about 1/3 of AML patients,
  the additional agent in treatment arm B was presumed to target this anomaly.
  All subjects had a FLT mutation, either internal tandem duplications (ITD)
  (divided into low vs high) +- mutations in the TKD domain, or TKD mutations
  only.  This was a stratification factor for treatment assignment in
  the study.   The levels of A, B, C correspond to increasing severity of the
  mutation burden.
}
\examples{
coxph(Surv(futime, death) ~ trt + flt3, data=myeloid)
# See the mstate vignette for a more complete analysis
}
\keyword{datasets}
\references{
  Le-Rademacher JG, Peterson RA, Therneau TM, Sanford BL, Stone RM,
  Mandrekar SJ.  Application of multi-state models in cancer clinical trials. 
  Clin Trials. 2018 Oct; 15 (5):489-498
  }