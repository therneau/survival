\name{cgd0}
\docType{data}
\alias{cgd0}
\title{Chronic Granulotomous Disease data}

\description{Data are from a placebo controlled trial of gamma
  interferon in chronic granulotomous disease (CGD).
  Contains the data on time to serious infections observed through
  end of study for each patient. 
}
\usage{cgd0}
\format{
  \describe{
   \item{id}{subject identification number}
   \item{center}{enrolling center }
   \item{random}{date of randomization }
   \item{treatment}{placebo or gamma interferon }
   \item{sex}{sex}
   \item{age}{age in years, at study entry }
   \item{height}{height in cm at study entry}
   \item{weight}{weight in kg at study entry}
   \item{inherit}{pattern of inheritance }
   \item{steroids}{use of steroids at study entry,1=yes}
   \item{propylac}{use of prophylactic antibiotics at study entry}
   \item{hos.cat}{a categorization of the centers into 4 groups}
   \item{futime}{days to last follow-up}
   \item{etime1-etime7}{up to 7 infection times for the subject}
   }
 }
 \details{
   This data set is in the form found in the references,
   with one line per patient and no recoding of the variables.
}
\note{
  The usage above assumes that the survival package is attached via
  \code{library(survival)} or via another package that depends on
  survival, e.g., mstate.  If not, then either \code{survival::cgd0}
  or \code{data(cgd, package="survival")} can be used.
}
\source{
  Fleming and Harrington, Counting Processes and Survival Analysis,
  appendix D.2. 
}
\seealso{\code{\link{cgd}}}
\keyword{datasets}
\keyword{survival}
