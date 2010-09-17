\name{rats}
\alias{rats}
\docType{data}
\title{Rat data from Gail et al.}
\description{48 rats were injected with a carcinogen, and then
  randomized to either drug or placebo.  The number of tumors ranges
  from 0 to 13; all rats were  censored at 6 months after randomization.
}
\usage{rats}
\format{
  \tabular{ll}{
    rat:\tab id\cr
    rx:\tab treatment,(1=drug, 0=control) \cr
    observation:\tab within rat\cr
    start:\tab entry time\cr
    stop:\tab exit time\cr
    status:\tab event status\cr
  }
}
\source{
  MH Gail, TJ Santner, and CC Brown (1980),
  An analysis of comparative carcinogenesis experiments based on
  multiple times to tumor.
  \emph{Biometrics} \bold{36}, 255--266.
}
\keyword{survival}
\keyword{datasets}
