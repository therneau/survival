\name{rats}
\alias{rats}
\docType{data}
\title{Rat treatment data from Mantel et al}
\description{Rat treatment data from Mantel et al.
  Three rats were chosen from each of 100 litters, one of which was
  treated with a drug, and then all followed for tumor incidence.
}
\usage{rats}
\format{
  \tabular{ll}{
    litter:\tab litter number from 1 to 100\cr
    rx:\tab treatment,(1=drug, 0=control) \cr
    time:\tab time to tumor or last follow-up\cr
    status:\tab event status, 1=tumor and 0=censored\cr
    sex:\tab male or female
  }
}
\source{
  N. R. Bohidar and J. L. Ciminera.
  Mantel-Haenszel analyses of litter-matched time to response data,
  with modifications for recovery of interlitter information.
  Cancer Research, 37:3863-3868, 1977. 5 
}
\references{
  E. W. Lee, L. J. Wei, and D. Amato,
    Cox-type regression analysis for large number of small groups of
    correlated failure time observations,
    in "Survival Analysis, State of the Art", Kluwer, 1992.
  }
\note{The subset of females (odd numbered litters) was used Lee et al.}
\keyword{survival}
\keyword{datasets}
