\name{ovarian}
\alias{ovarian}
\docType{data}
\title{Ovarian Cancer Survival Data}
\usage{ovarian
data(cancer, package="survival")
}
\description{Survival in a randomised trial comparing two treatments for
    ovarian cancer}
  \format{
    A data frame with 26 observations and the following variables
    \describe{
      \item{\code{futime}}{survival or censoring time}
      \item{\code{fustat}}{0=censored, 1=death}
      \item{\code{age}}{age in years}
      \item{\code{resid.ds}}{residual disease present (1=no,2=yes)}
      \item{\code{rx}}{treatment group, 1= cyclophosphamide,
	2= cyclophosphamide + adriamycin}
      \item{\code{ecog.ps}}{ECOG performance status (1 is better)}
    }
  }
}
\references{
   Edmunson, J.H., Fleming, T.R., Decker, D.G.,
    Malkasian, G.D., Jefferies, J.A., Webb, M.J., and Kvols, L.K.,
    Different Chemotherapeutic Sensitivities and Host Factors Affecting
    Prognosis in Advanced Ovarian Carcinoma vs. Minimal Residual Disease.
    Cancer Treatment Reports, 63:241-47, 1979.
  }
\keyword{datasets}
\keyword{survival}
