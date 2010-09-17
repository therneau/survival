\name{veteran}
\alias{veteran}
\docType{data}
\title{Veterans' Administration Lung Cancer study}
\description{Randomised trial of two treatment regimens for lung cancer.
  This is a standard survival analysis data set.} 
\usage{veteran}
\format{
  \tabular{ll}{
    trt:\tab 1=standard 2=test\cr
    celltype:\tab 1=squamous,  2=smallcell,  3=adeno,  4=large\cr
    time:\tab survival time\cr
    status:\tab censoring status\cr
    karno:\tab Karnofsky performance score (100=good)\cr
    diagtime:\tab months from diagnosis to randomisation\cr
    age:\tab in years\cr
    prior:\tab prior therapy 0=no, 1=yes\cr
  }
}
\source{
  D Kalbfleisch and RL Prentice (1980),
  \emph{The Statistical Analysis of Failure Time Data}.
  Wiley, New York.
}
\keyword{datasets}
