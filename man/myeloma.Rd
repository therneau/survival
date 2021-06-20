\name{myeloma}
\alias{myeloma}
\docType{data}
\title{
Survival times of patients with multiple myeloma
}
\description{
Survival times of 3882 subjects with multiple myeloma, seen at Mayo
Clinic from 1947--1996.
}
\usage{myeloma
data("cancer", package="survival")}
\format{
  A data frame with 3882 observations on the following 5 variables.
  \describe{
    \item{\code{id}}{subject identifier}
    \item{\code{year}}{year of entry into the study}
    \item{\code{entry}}{time from diagnosis of MM until entry (days)}
    \item{\code{futime}}{follow up time (days)}
    \item{\code{death}}{status at last follow-up: 0 = alive, 1 = death}
  }
}
\details{
Subjects who were diagnosed at Mayo will have \code{entry} =0, those who
were diagnosed elsewhere and later referred will have positive values.
}
\references{
  R. Kyle, Long term survival in multiple myeloma.
  New Eng J Medicine, 1997
}
\examples{
# Incorrect survival curve, which ignores left truncation
fit1 <- survfit(Surv(futime, death) ~ 1, myeloma)
# Correct curve
fit2 <- survfit(Surv(entry, futime, death) ~1, myeloma)
}
\keyword{datasets}
