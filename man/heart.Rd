\name{heart}
\docType{data}
\alias{jasa1}
\alias{jasa}
\alias{heart}
\title{Stanford Heart Transplant data}
\description{Survival of patients on the waiting list for the Stanford
  heart transplant program.}
\usage{
heart
jasa
jasa1
}
\format{
  jasa: original data
  \tabular{ll}{
   birth.dt:\tab birth date \cr
   accept.dt:\tab acceptance into program \cr
   tx.date:\tab transplant date \cr
   fu.date:\tab end of followup \cr
   fustat:\tab dead or alive \cr
   surgery:\tab prior bypass surgery\cr
   age: \tab age (in days)\cr
   futime:\tab followup time\cr
   wait.time:\tab time before transplant\cr
   transplant:\tab transplant indicator\cr
   mismatch:\tab mismatch score\cr
   hla.a2:\tab particular type of mismatch\cr
   mscore:\tab another mismatch score\cr
   reject:\tab rejection occurred\cr
 }

 jasa1, heart: processed data
  \tabular{ll}{
    start, stop, event: \tab Entry and exit time and status for this interval of time\cr
    age:\tab age-48 years\cr
    year:\tab year of acceptance (in years after 1 Nov 1967)\cr
    surgery:\tab prior bypass surgery 1=yes\cr
    transplant: \tab received transplant 1=yes\cr
    id:\tab patient id\cr
  }
}
\seealso{\code{\link{stanford2}}}
\source{
  J Crowley and M Hu (1977),
  Covariance analysis of heart transplant survival data.
  \emph{Journal of the American Statistical Association},
  \bold{72}, 27--36.
}
\keyword{datasets}
\keyword{survival}
