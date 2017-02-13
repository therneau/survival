\name{pbc}
\alias{pbc}
\docType{data}
\title{Mayo Clinic Primary Biliary Cirrhosis Data}
\description{D
This data is from the Mayo Clinic
trial in primary biliary cirrhosis (PBC) of the liver conducted between 1974
and 1984.  
A total of 424 PBC patients, referred to Mayo Clinic during that ten-year
interval, met eligibility criteria for the randomized placebo controlled
trial of the drug D-penicillamine.  The first 312 cases in the data set
participated in the randomized trial and contain largely complete data.  The
additional 112 cases did not participate in the clinical trial, but consented
to have basic measurements recorded and to be followed for survival.  Six of
those cases were lost to follow-up shortly after diagnosis, so the data here
are on an additional 106 cases as well as the 312 randomized participants.
    
A nearly identical data set found in appendix D of Fleming and Harrington;
this version has fewer missing values.
}
\usage{pbc}
\format{
  \tabular{ll}{
    age:\tab in years\cr
    albumin:\tab serum albumin (g/dl)\cr
    alk.phos:\tab alkaline phosphotase (U/liter)\cr
    ascites:\tab presence of ascites \cr
    ast:\tab aspartate aminotransferase, once called SGOT (U/ml)\cr
    bili:\tab serum bilirunbin (mg/dl)\cr
    chol:\tab serum cholesterol (mg/dl)\cr
    copper:\tab urine copper (ug/day)\cr
    edema:\tab 0 no edema, 0.5 untreated or successfully treated\cr
    \tab 1 edema despite diuretic therapy\cr
    hepato:\tab presence of hepatomegaly or enlarged liver\cr
    id:\tab case number\cr
    platelet:\tab platelet count\cr
    protime:\tab standardised blood clotting time\cr
    sex:\tab m/f\cr
    spiders:\tab blood vessel malformations in the skin\cr
    stage:\tab histologic stage of disease (needs biopsy)\cr
    status:\tab status at endpoint, 0/1/2 for censored, transplant, dead\cr
    time: \tab number of days between registration and the earlier of death,\cr
       \tab transplantion, or study analysis in July, 1986\cr
    trt:\tab 1/2/NA for D-penicillmain, placebo, not randomised\cr
    trig:\tab triglycerides (mg/dl)\cr
  }
}
\source{
  T Therneau and P Grambsch (2000),
  \emph{Modeling Survival Data: Extending the Cox Model},
  Springer-Verlag, New York.
  ISBN: 0-387-98784-3.
}
\keyword{datasets}
