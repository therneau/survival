\name{udca}
\alias{udca}
\alias{udca1}
\alias{udca2}
\docType{data}
\title{Data from a trial of usrodeoxycholic acid
}
\description{
  Data from a trial of ursodeoxycholic acid (UDCA) in patients with primary
  biliary cirrohosis (PBC).  
}
\usage{udca
udca2
data(udca, package="survival")
}
\format{
  A data frame with 170 observations on the following 15 variables.
  \describe{
    \item{\code{id}}{subject identifier}
    \item{\code{trt}}{treatment of 0=placebo, 1=UDCA}
    \item{\code{entry.dt}}{date of entry into the study}
    \item{\code{last.dt}}{date of last on-study visit}
    \item{\code{stage}}{stage of disease}
    \item{\code{bili}}{bilirubin value at entry}
    \item{\code{riskscore}}{the Mayo PBC risk score at entry}
    \item{\code{death.dt}}{date of death}
    \item{\code{tx.dt}}{date of liver transplant}
    \item{\code{hprogress.dt}}{date of histologic progression}
    \item{\code{varices.dt}}{appearance of esphogeal varices}
    \item{\code{ascites.dt}}{appearance of ascites}
    \item{\code{enceph.dt}}{appearance of encephalopathy}
    \item{\code{double.dt}}{doubling of initial bilirubin}
    \item{\code{worsen.dt}}{worsening of symptoms by two stages}
  }
}
\details{
  This data set is used in the Therneau and Grambsh.  The \code{udca1}
  data set contains the baseline variables along with the time until the
  first endpoint (any of death, transplant, \ldots, worsening).
  The \code{udca2} data set treats all of the endpoints as parallel
  events and has a stratum for each.
}

\references{
T. M. Therneau and P. M. Grambsch, Modeling survival data: extending the Cox
model.  Springer, 2000.

K. D. Lindor, E. R. Dickson, W. P Baldus, R.A. Jorgensen, J. Ludwig,
P. A. Murtaugh, J. M. Harrison, R. H. Weisner, M. L. Anderson,
S. M. Lange, G. LeSage, S. S. Rossi and A. F. Hofman.
Ursodeoxycholic acid in the treatment of primary biliary cirrhosis.
Gastroenterology, 106:1284-1290, 1994.
}
\examples{
# values found in table 8.3 of the book
fit1 <- coxph(Surv(futime, status) ~ trt + log(bili) + stage,
          cluster =id , data=udca1)
fit2 <- coxph(Surv(futime, status) ~ trt + log(bili) + stage +
          strata(endpoint), cluster=id,  data=udca2)

}
\keyword{datasets}
