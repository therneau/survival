\name{stanford2}
\alias{stanford2}
\docType{data}
\title{More Stanford Heart Transplant data}
\description{
  This contains the Stanford Heart Transplant data in a different
  format.  The main data set is in \code{\link{heart}}.
}
\usage{stanford2}
\format{
  \tabular{ll}{
    id:  \tab ID number\cr
    time:\tab survival or censoring time\cr
    status:\tab censoring status\cr
    age: \tab in years\cr
    t5: \tab T5 mismatch score\cr
  }
}
\seealso{
  \code{\link{predict.survreg}},
  \code{\link{heart}}
}
\source{
  LA Escobar and WQ Meeker Jr (1992),
  Assessing influence in regression analysis with censored data.
  \emph{Biometrics} \bold{48}, 507--528.
  Page 519.
}
\keyword{datasets}
\keyword{survival}
