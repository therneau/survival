\name{cgd}
\docType{data}
\alias{cgd}
\alias{cgd.raw}
\title{Chronic Granulotamous Disease data}

\description{Data are from a placebo controlled trial of gamma
  interferon in chronic granulotomous disease (CGD).
  Contains the data on time to serious infections observed through
  end of study for each patient. 
}
\usage{cgd}
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
   \item{tstart, tstop}{start and end of each time interval }
   \item{status}{1=the interval ends with an infection }
   \item{enum}{observation number within subject}
   }
 }
 \details{
   The \code{cgd0} data set is in the form found in the references,
   with one line per patient and no recoding of the variables.
   The \code{cgd} data set (this one) has been cast into (start, stop]
   format with one line per event, and covariates
   such as center recoded as factors
   to include meaningful labels. 
}
\source{
  Fleming and Harrington, Counting Processes and Survival Analysis,
  appendix D.2. 
}
\seealso{\code{link{cgd0}}}
\keyword{datasets}
\keyword{survival}
