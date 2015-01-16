\name{cgd}
\docType{data}
\alias{cgd}
\alias{cgd.raw}
\title{Chronic Granulotomous Disease data}

\description{Data are from a placebo controlled trial of gamma
interferon in chronic granulotomous disease (CGD). Uses the complete
data on time to first serious infection observed through end of study
for each patient, which includes the initial serious infections observed
through the 7/15/89 interim analysis data cutoff, plus the residual data
on occurence of initial serious infections between the interim analysis
cutoff and the final blinded study visit for each patient. Only one
patient was taken off on the day of his last infection.
}
\usage{cgd}
\format{
  \describe{
   \item{id}{subject identifiction number}
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
\usage{cgdraw}
\format{
  \describe{
   \item{id}{subject identifiction number}
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
   The \code{cgdraw} data set is in the form found in the references,
   with one line per patient and no recoding of the variables.
   The \item{cgd} data set has been further processed so as to have one
   line per event, with covariates such as center recoded as factors
   to include meaningful labels. 

   The second data set is ready for use in \code{coxph} calls, and will
   be found in examples.
   The raw data set is used in one vignette to show how to create
   multi-event data sets.
}
\source{
  Fleming and Harrington, Counting Processes and Survival Analysis,
  appendix D.2. 
}
\keyword{datasets}
\keyword{survival}
