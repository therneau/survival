\name{colon}
\alias{colon}
\title{Chemotherapy for Stage B/C colon cancer}
\usage{colon}
\description{These are data from one of the first successful trials of
    adjuvant chemotherapy for colon cancer. Levamisole is a low-toxicity
    compound previously used to treat worm infestations in animals; 5-FU
    is a moderately toxic (as these things go) chemotherapy agent. There
  are two records per person, one for recurrence and one for death}
\format{
\tabular{ll}{
id:\tab id\cr
study:\tab 1 for all patients\cr
rx:\tab Treatment - Obs(ervation), Lev(amisole), Lev(amisole)+5-FU\cr
sex:\tab  1=male\cr
age:\tab in years\cr
obstruct:\tab obstruction of colon by tumour\cr
perfor:\tab perforation of colon\cr
adhere:\tab adherence to nearby organs\cr
nodes:\tab number of lymph nodes with detectable cancer\cr
time:\tab days until event or censoring\cr
status:\tab censoring status\cr
differ:\tab differentiation of tumour (1=well, 2=moderate, 3=poor)\cr
extent:\tab Extent of local spread (1=submucosa, 2=muscle, 3=serosa,
4=contiguous structures)\cr
surg:\tab time from surgery to registration (0=short, 1=long)\cr
node4:\tab more than 4 positive lymph nodes\cr
etype:\tab event type: 1=recurrence,2=death\cr
}}

\source{Danyu Lin}
\keyword{survival}
