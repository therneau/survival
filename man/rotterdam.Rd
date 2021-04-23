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
These data sets are used in the paper by Royston and Altman.
The Rotterdam data is used to create a fitted model, and the GBSG data for 
validation of the model.  The paper gives references for the data source.
}
\seealso{
  \code{\link{gbsg}}
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
