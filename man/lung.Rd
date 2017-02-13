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
\usage{
lung
cancer
}
\format{
  \tabular{ll}{
    inst:\tab Institution code\cr
    time:\tab Survival time in days\cr
    status:\tab censoring status 1=censored, 2=dead\cr
    age:\tab Age in years\cr
    sex:\tab  Male=1 Female=2\cr
    ph.ecog:\tab ECOG performance score (0=good 5=dead)\cr
    ph.karno:\tab Karnofsky performance score (bad=0-good=100) rated by physician\cr
    pat.karno:\tab Karnofsky performance score as rated by patient\cr
    meal.cal:\tab Calories consumed at meals\cr
    wt.loss:\tab Weight loss in last six months\cr
  }
}
\source{Terry Therneau}
\references{
  Loprinzi CL. Laurie JA. Wieand HS. Krook JE. Novotny PJ.
  Kugler JW. Bartel J. Law M. Bateman M. Klatt NE. et al.
  Prospective evaluation of prognostic variables from patient-completed
  questionnaires. North Central Cancer Treatment Group.
  Journal of Clinical Oncology. 12(3):601-7, 1994. }
\keyword{datasets}
