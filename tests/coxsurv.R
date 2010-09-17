options(na.action=na.exclude) # preserve missings
options(contrasts=c('contr.treatment', 'contr.poly')) #ensure constrast type
library(survival)

#
# Test out subscripting in the case of a coxph survival curve
#
aeq <- function(x,y, ...) all.equal(as.vector(x), as.vector(y), ...)

fit <- coxph(Surv(time, status) ~ age + sex + meal.cal + strata(ph.ecog),
		data=cancer)

surv1 <- survfit(fit)
temp <- surv1[2:3]

which <- cumsum(surv1$strata)
zed   <- (which[1]+1):(which[3])
aeq(surv1$surv[zed], temp$surv)
aeq(surv1$time[zed], temp$time)

#
# Now a result with a matrix of survival curves
#
dummy <- data.frame(age=c(30,40,60), sex=c(1,2,2), meal.cal=c(500, 1000, 1500))
surv2 <- survfit(fit, newdata=dummy)

zed <- 1:which[1]
aeq(surv2$surv[zed,1], surv2[1,1]$surv)
aeq(surv2$surv[zed,2], surv2[1,2]$surv)
aeq(surv2$surv[zed,3], surv2[1,3]$surv)
aeq(surv2$surv[zed, ], surv2[1,1:3]$surv)
aeq(surv2$surv[zed],   (surv2[1]$surv)[,1])
aeq(surv2$surv[zed, ], surv2[1, ]$surv)

rm(fit, surv1, temp, which, zed,dummy, surv2)
