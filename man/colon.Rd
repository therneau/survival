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
\note{The study is originally described in Laurie (1989). 
  The main report is found in Moertel (1990).  This data set is closest
  to that of the final report in Moertel (1991).
   A version of the data with less follow-up time was used in the
   paper by Lin (1994).
   }
\references{
  JA Laurie, CG Moertel, TR Fleming, HS Wieand, JE Leigh, J Rubin,
  GW McCormack, JB Gerstner, JE Krook and J Malliard.  Surgical
  adjuvant therapy of large-bowel carcinoma: An evaluation of
  levamisole and the combination of levamisole and fluorouracil:
  The North Central Cancer Treatment Group and the Mayo Clinic.
  J Clinical Oncology, 7:1447-1456, 1989.

DY Lin.  Cox regression analysis of multivariate failure time data:
  the marginal approach.  Statistics in Medicine, 13:2233-2247, 1994.

  CG Moertel, TR Fleming, JS MacDonald, 
DG Haller, JA Laurie, PJ Goodman, JS Ungerleider, 
WA Emerson, DC Tormey, JH Glick, MH Veeder and JA Maillard. 
Levamisole and fluorouracil for adjuvant therapy of 
resected colon carcinoma. New England J of Medicine, 332:352-358,
1990.

 CG Moertel,  TR Fleming,  JS MacDonald, DG Haller, JA Laurie, CM Tangen,
 JS Ungerleider, WA Emerson, DC Tormey, JH Glick, MH Veeder and JA
 Maillard, Fluorouracil plus Levamisole as and effective adjuvant
 therapy after resection of stage II colon carcinoma: a final report.
 Annals of Internal Med, 122:321-326, 1991.
  }
\keyword{survival}
