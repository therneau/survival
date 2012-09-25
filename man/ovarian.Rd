\name{ovarian}
\alias{ovarian}
\docType{data}
\title{Ovarian Cancer Survival Data}
\usage{ovarian}
\description{Survival in a randomised trial comparing two treatments for
    ovarian cancer}
\format{
  \tabular{ll}{
    futime:\tab survival or censoring time\cr
    fustat:\tab censoring status\cr
    age: \tab in years\cr
    resid.ds:\tab residual disease present (1=no,2=yes)\cr
    rx:\tab treatment group\cr
    ecog.ps:\tab ECOG performance status (1 is better, see reference)\cr
  }
}
\source{Terry Therneau}
\references{
   Edmunson, J.H., Fleming, T.R., Decker, D.G.,
    Malkasian, G.D., Jefferies, J.A., Webb, M.J., and Kvols, L.K.,
    Different Chemotherapeutic Sensitivities and Host Factors Affecting
    Prognosis in Advanced Ovarian Carcinoma vs. Minimal Residual Disease.
    Cancer Treatment Reports, 63:241-47, 1979.
    
  for ECOG performance status:
  \url{http://ecog.org/general/perf_stat.html}
  }
\keyword{datasets}
\keyword{survival}
