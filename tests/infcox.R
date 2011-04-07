options(na.action=na.exclude) # preserve missings
options(contrasts=c('contr.treatment', 'contr.poly')) #ensure constrast type
library(survival)

#
# A test to exercise the "infinity" check on 2 variables
#
test3 <- data.frame(futime=1:12, fustat=c(1,0,1,0,1,0,0,0,0,0,0,0),
		   x1=rep(0:1,6), x2=c(rep(0,6), rep(1,6)))
coxph(Surv(futime, fustat) ~ x1 + x2, test3, iter=25)
