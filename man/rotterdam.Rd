\name{rotterdam}
\alias{rotterdam}
\docType{data}
\title{Breast cancer data set used in Royston and Altman (2013)}
\description{
The \code{rotterdam} data set includes 2982 primary breast cancers patients
whose data whose records were included in the Rotterdam tumor bank.
}
\usage{rotterdam
data(cancer, package="survival")
}
\format{
  A data frame with 2982 observations on the following 15 variables.
  \describe{
    \item{\code{pid}}{patient identifier}
    \item{\code{year}}{year of surgery}
    \item{\code{age}}{age at surgery}
    \item{\code{meno}}{menopausal status (0= premenopausal, 1= postmenopausal)}
    \item{\code{size}}{tumor size, a factor with levels \code{<=20} \code{20-50} \code{>50}}
    \item{\code{grade}}{differentiation grade}
    \item{\code{nodes}}{number of positive lymph nodes}
    \item{\code{pgr}}{progesterone receptors (fmol/l)}
    \item{\code{er}}{estrogen receptors (fmol/l)}
    \item{\code{hormon}}{hormonal treatment (0=no, 1=yes)}
    \item{\code{chemo}}{chemotherapy}
    \item{\code{rtime}}{days to relapse or last follow-up}
    \item{\code{recur}}{0= no relapse, 1= relapse}
    \item{\code{dtime}}{days to death or last follow-up}
    \item{\code{death}}{0= alive, 1= dead}
  }
}
\details{
  These data sets are used in the paper by Royston and Altman that is
  referenced below.
The Rotterdam data is used to create a fitted model, and the GBSG data for 
validation of the model.  The paper gives references for the data
source.

There are 43 subjects who have died without recurrence, but whose death
time is greater than the censoring time for recurrence.
A common way that this happens is that a death date is updated in the
health record sometime after the research study ended, and said value
is then picked up when a study data set is created.
But it raises serious questions about censoring.
For instance subject 40 is censored for recurrence at 4.2 years and died
at 6.6 years; when creating the endpoint of recurrence free survival
(earlier of recurrence or death), treating them as a death at 6.6 years
implicitly assumes that they were recurrence free just before death.
For this to be true we would have to assume that if they had progressed in
the 2.4 year interval before death (off study),
that this information would also have been noted
in their general medical record, and would also be captured in
the study data set.
However, that may be unlikely.  Death information is often in a
centralized location in electronic health records, easily accessed by a
programmer and merged with the study data, while recurrence may
require manual review.  How best to address this is an open issue.
}
\seealso{
  \code{\link{gbsg}}
}
\references{
Patrick Royston and Douglas Altman, External validation of a Cox prognostic
model: principles and methods.  BMC Medical Research Methodology 2013, 13:33
}
\examples{
status  <- pmax(rotterdam$recur, rotterdam$death)
rfstime <- with(rotterdam, ifelse(recur==1, rtime, dtime))
fit1 <- coxph(Surv(rfstime, status) ~ pspline(age) + meno + size + 
        pspline(nodes) + er,
     data=rotterdam, subset = (nodes > 0))
# Royston and Altman used fractional polynomials for the nonlinear terms
}
\keyword{datasets}
\keyword{survival}
