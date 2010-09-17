\name{cgd}
\docType{data}
\alias{cgd}
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
  \tabular{ll}{
   id:\tab subject identifiction number\cr
   center:\tab enrolling center \cr
   random:\tab date of randomization \cr
   treatment:\tab placebo or gamma interferon \cr
   sex:\tab  \cr
   age:\tab age in years, at study entry \cr
   height:\tab height in cm at study entry\cr
   weight: \tab weight in kg at study entry\cr
   inherit:\tab pattern of inheritance \cr
   steroids:\tab use of steroids at study entry,1=yes\cr
   propylac: \tab use of prophylactic antibiotics at study entry\cr
   hos.cat:\tab a categorization of the centers into 4 groups\cr
   tstart, tstop:\tab start and end of each time interval \cr
   status:\tab 1=the interval ends with an infection \cr
   enum: \tab observation number within subject\cr
   }
 }
\source{
  Fleming and Harrington, Counting Processes and Survival Analysis,
  appendix D.2. 
}
\keyword{datasets}
\keyword{survival}
