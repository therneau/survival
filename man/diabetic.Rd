\name{diabetic}
\alias{diabetic}
\docType{data}
\title{Ddiabetic retinopathy}
\description{
  Partial results from a trial of laser coagulation for the treatment
  of diabetic retinopathy.
}
\usage{diabetic}
\format{
  A data frame with 394 observations on the following 8 variables.
  \describe{
    \item{\code{id}}{subject id}
    \item{\code{laser}}{laser type: \code{xenon} or \code{argon}}
    \item{\code{age}}{age at diagnosis}
    \item{\code{eye}}{a factor with levels of \code{left} \code{right}}
    \item{\code{trt}}{treatment: 0 = no treatment, 1= laser}
    \item{\code{risk}}{risk group of 6-12}
    \item{\code{time}}{time to event or last follow-up}
    \item{\code{status}}{status of 0= censored or 1 = visual loss}
  }
}
\details{
The 197 patients in this dataset were a 50\% random sample of the
patients with "high-risk" diabetic retinopathy as defined by the
Diabetic Retinopathy Study (DRS).  Each patient had one eye randomized
to laser treatment and the other eye received no treatment.  For each
eye, the event of interest was the time from initiation of treatment
to the time when visual acuity dropped below 5/200 two visits in a row.
Thus there is a built-in lag time of
approximately 6 months (visits were every 3 months).  Survival times
in this dataset are therefore the actual time to blindness in months,
minus the minimum possible time to event (6.5 months).  Censoring was
caused by death, dropout, or end of the study.
}
\references{
  Huster, Brookmeyer and Self, Biometrics, 1989.

  American Journal of Ophthalmology, 1976, 81:4, pp 383-396
}
\examples{
# juvenile diabetes is defined as and age less than 20
juvenile <- 1*(diabetic$age < 20)
coxph(Surv(time, status) ~ trt + juvenile, cluster= id,
            data= diabetic)
}
\keyword{datasets}
\keyword{survival}