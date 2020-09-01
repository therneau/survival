\name{breastca}
\alias{rotterdam}
\alias{gbsg}
\docType{data}
\title{Breast cancer data sets used in Royston and Altman (2013)}
\description{
The \code{rotterdam} data set includes 2982 primary breast cancers patients
whose data whose records were included in the Rotterdam tumor bank.
The \code{gbsg} data set contains patient records from a 1984-1989 trial 
conducted by the German Breast Cancer Study Group (GBSG) of 720 patients
with node positive breast cancer; it retains the 686 patients with 
complete data for the prognostic variables.
}
\usage{data("rotterdam")}
\format{
  A data frame with 2982 observations on the following 15 variables.
  \describe{
    \item{\code{pid}}{patient identifier}
    \item{\code{year}}{year of cancer incidence}
    \item{\code{age}}{age}
    \item{\code{meno}}{menopausal status (0= premenopausal, 1= postmenopausal)}
    \item{\code{size}}{tumor size, a factor with levels \code{<=20} \code{20-50} \code{>50}}
    \item{\code{grade}}{tumor grade}
    \item{\code{nodes}}{number of positive lymph nodes}
    \item{\code{pgr}}{progesterone receptors (fmol/l)}
    \item{\code{er}}{estrogen receptors (fmol/l)}
    \item{\code{hormon}}{hormonal treatment (0=no, 1=yes)}
    \item{\code{chemo}}{chemotherapy}
    \item{\code{rtime}}{days to recurrence or last follow-up}
    \item{\code{recur}}{0= no recurrence, 1= recurrence}
    \item{\code{dtime}}{days to death or last follow-up}
    \item{\code{death}}{0= alive, 1= dead}
  }
}
\usage{data("gbsg")}
\format{
   A data set with 686 observations and 11 variables.
   \describe{
   \item{\code{pid}}{patient identifier}
    \item{\code{age}}{age}
    \item{\code{meno}}{menopausal status (0= premenopausal, 1= postmenopausal)}
    \item{\code{size}}{tumor size}
    \item{\code{grade}}{tumor grade}
    \item{\code{nodes}}{number of positive lymph nodes}
    \item{\code{pgr}}{progesterone receptors (fmol/l)}
    \item{\code{er}}{estrogen receptors (fmol/l)}
    \item{\code{hormon}}{hormonal therapy, 0= no, 1= yes}
    \item{\code{rfstime}}{recurrnece free survival time; days to first of reccurence, death or last follow-up}
    \item{\code{status}}{0= alive without recurrence, 1= recurrence or death}
\details{
These data sets are used in the paper by Royston and Altman.
The Rotterdam data is used to create a fitted model, and the GBSG data for 
validation of the model.  The paper gives references for the data source.
}
\references{
Patrick Royston and Douglas Altman, External validation of a Cox prognostic
model: principles and methods.  BMC Medical Research Methodology 2013, 13:33
}
\examples{
rfstime <- pmin(rotterdam$rtime, rotterdam$dtime)
status  <- pmax(rotterdam$recur, rotterdam$death)
fit1 <- coxph(Surv(rfstime, status) ~ pspline(age) + meno + size + 
        pspline(nodes) + er,
     data=rotterdam, subset = (nodes > 0))
# Royston and Altman used fractional polynomials for the nonlinear terms
}
\keyword{datasets}
\keyword{survival}
