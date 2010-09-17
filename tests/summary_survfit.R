## check that the scale option to summary.survfit works
##  Marc Schwartz reported this as a bug in 2.35-3.
library(survival)
summary( survfit( Surv(futime, fustat)~1, data=ovarian))
summary( survfit( Surv(futime, fustat)~1, data=ovarian), scale=365.25)
