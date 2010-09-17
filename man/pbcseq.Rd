\name{pbcseq}
\alias{pbcseq}
\docType{data}
\title{Mayo Clinic Primary Biliary Cirrhosis, sequential data}
\description{
This data is a continuation of the PBC data set,
and contains the follow-up laboratory data for each study patient.
An analysis based on the data can be found in Murtagh, et. al.
  
The primary PBC data set contains only baseline measurements 
of the laboratory 
paramters.  This data set contains multiple laboratory results, but
only on the 312 randomized patients.  Some baseline data values in this file
differ from the original PBC file, for instance, the data errors in
prothrombin time and age which were discovered after the orignal analysis
(see Fleming and Harrington, figure 4.6.7).
   
One "feature" of the data deserves special comment.  The last
observation before death or liver transplant often has many more
missing covariates than other data rows.  The original clinical
protocol for these patients specified visits at 6 months, 1 year, and
annually thereafter.  At these protocol visits lab values were
obtained for a large pre-specified battery of tests.  "Extra" visits,
often undertaken because of worsening medical condition, did not
necessarily have all this lab work. 
The missing values are thus potentially informative.
}
\usage{pbc}
\format{
  \tabular{ll}{
    id:\tab case number\cr
    age:\tab in years\cr
    sex:\tab m/f\cr
    trt:\tab 1/2/NA for D-penicillmain, placebo, not randomised\cr
    time:\tab number of days between registration and the earlier of death,\cr
       \tab transplantion, or study analysis in July, 1986\cr
    status:\tab status at endpoint, 0/1/2 for censored, transplant, dead\cr
    day:\tab number of days between enrollment and this visit date\cr
       \tab all measurements below refer to this date\cr
    albumin:\tab serum albumin (mg/dl)\cr
    alk.phos:\tab alkaline phosphotase (U/liter)\cr
    ascites:\tab presence of ascites \cr
    ast:\tab aspartate aminotransferase, once called SGOT (U/ml)\cr
    bili:\tab serum bilirunbin (mg/dl)\cr
    chol:\tab serum cholesterol (mg/dl)\cr
    copper:\tab urine copper (ug/day)\cr
    edema:\tab 0 no edema, 0.5 untreated or successfully treated\cr
    \tab 1 edema despite diuretic therapy\cr
    hepato:\tab presence of hepatomegaly or enlarged liver\cr
    platelet:\tab platelet count\cr
    protime:\tab standardised blood clotting time\cr
    spiders:\tab blood vessel malformations in the skin\cr
    stage:\tab histologic stage of disease (needs biopsy)\cr
    trig:\tab triglycerides (mg/dl)\cr
  }
}
\source{
  T Therneau and P Grambsch,
  "Modeling Survival Data: Extending the Cox Model",
  Springer-Verlag, New York, 2000.
  ISBN: 0-387-98784-3.
}
\references{
  Murtaugh PA. Dickson ER. Van Dam GM. Malinchoc M. Grambsch PM.
  Langworthy AL. Gips CH.  "Primary biliary cirrhosis: prediction of short-term
  survival based on repeated patient visits." Hepatology. 20(1.1):126-34, 1994.

  Fleming T and Harrington D., "Counting Processes and Survival Analysis",
  Wiley, New York, 1991.
  }
\keyword{datasets}
