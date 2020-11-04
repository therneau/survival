\name{lung}
\docType{data}
\alias{cancer}
\alias{lung}
\title{NCCTG Lung Cancer Data}
\description{
  Survival in patients with advanced lung cancer from the North
  Central Cancer Treatment Group.  Performance
  scores rate how well the patient can perform usual daily activities.
}
\usage{lung}
\format{
  \tabular{ll}{
    inst:\tab Institution code\cr
    time:\tab Survival time in days\cr
    status:\tab censoring status 1=censored, 2=dead\cr
    age:\tab Age in years\cr
    sex:\tab  Male=1 Female=2\cr
    ph.ecog:\tab ECOG performance score as rated by the physician.
    0=asymptomatic, 1= symptomatic but completely ambulatory, 2= in bed
    <50\% of the day, 3= in bed > 50\% of the day but not bedbound, 4 =
    bedbound\cr
    ph.karno:\tab Karnofsky performance score (bad=0-good=100) rated by physician\cr
    pat.karno:\tab Karnofsky performance score as rated by patient\cr
    meal.cal:\tab Calories consumed at meals\cr
    wt.loss:\tab Weight loss in last six months\cr
  }
}
\note{
  The use of 1/2 for alive/dead instead of the usual 0/1 is a historical
  footnote.
  For data contained on punch cards, IBM 360 Fortran treated blank as a zero,
  which led to a policy within the section of Biostatistics to never
  use "0" as a data value since one could not distinguish it from a
  missing value.
  The policy became a habit, as is often the case; and the 1/2 coding
  endured long beyond the demise of punch cards and Fortran.
}
  

\source{Terry Therneau}
\references{
  Loprinzi CL. Laurie JA. Wieand HS. Krook JE. Novotny PJ.
  Kugler JW. Bartel J. Law M. Bateman M. Klatt NE. et al.
  Prospective evaluation of prognostic variables from patient-completed
  questionnaires. North Central Cancer Treatment Group.
  Journal of Clinical Oncology. 12(3):601-7, 1994. }
\keyword{datasets}
