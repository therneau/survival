options(na.action=na.exclude) # preserve missings
options(contrasts=c('contr.treatment', 'contr.poly')) #ensure constrast type
library(survival)

#
# A test to exercise the "infinity" check on 2 variables
#
test3 <- data.frame(futime=1:12, fustat=c(1,0,1,0,1,0,0,0,0,0,0,0),
		   x1=rep(0:1,6), x2=c(rep(0,6), rep(1,6)))

# This will produce a warning message, which is the point of the test.
# The variance is close to singular and gives different answers 
#  on different machines
fit3 <- coxph(Surv(futime, fustat) ~ x1 + x2, test3, iter=25)

all.equal(round(fit3$coef,1), c(-22.2, -22.9), check.attributes=FALSE)
all.equal(round(fit3$log, 4),c(-6.8669, -1.7918))
