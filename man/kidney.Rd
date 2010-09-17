\name{kidney}
\alias{kidney}
\title{Kidney catheter data}
\format{
\tabular{ll}{
patient:\tab id\cr
time:\tab time\cr
status:\tab event status\cr
age:\tab in years\cr
sex:\tab 1=male, 2=female\cr
disease:\tab  disease type (0=GN, 1=AN, 2=PKD, 3=Other)\cr
frail:\tab frailty estimate from original paper\cr
}}
\description{
Data on the recurrence times to infection, at the point of insertion of
  the catheter, for kidney patients using portable dialysis equipment.
  Catheters may be removed for reasons other than infection, in which case
  the observation is censored.  Each patient has exactly 2 observations.

This data has often been used to illustrate the use of random effects
(frailty) in a survival model.  However, one of the males (id 21) is a
large outlier, with much longer survival than his peers.  If this
observation is removed no evidence remains for a random subject effect.
  
}
\section{Note}{
  The original paper ignored the issue of tied times and so is not
  exactly reproduced by the survival package.
}
\examples{
kfit <- coxph(Surv(time, status)~ age + sex + disease + frailty(id), kidney)
kfit0 <- coxph(Surv(time, status)~ age + sex + disease, kidney)
kfitm1 <- coxph(Surv(time,status) ~ age + sex + disease + 
		frailty(id, dist='gauss'), kidney)

}
\source{
  CA McGilchrist, CW Aisbett (1991),
  Regression with frailty in survival analysis. 
  \emph{Biometrics} \bold{47}, 461--66.
}
\keyword{survival}
