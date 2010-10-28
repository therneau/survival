options(na.action=na.exclude) # preserve missings
options(contrasts=c('contr.treatment', 'contr.poly')) #ensure constrast type
library(survival)

#
# Test out the revised model.matrix code
#
test1 <- data.frame(time=  c(9, 3,1,1,6,6,8),
                    status=c(1,NA,1,0,1,1,0),
                    x=     c(0, 2,1,1,1,0,0),
                    z= factor(c('a', 'a', 'b', 'b', 'c', 'c', 'a')))

fit1 <- coxph(Surv(time, status) ~ z, test1, iter=1)
fit2 <- coxph(Surv(time, status) ~z, test1, x=T, iter=1)
all.equal(model.matrix(fit1), fit2$x)

# This has no level 'b', make sure dummies recode properly
test2 <- data.frame(time=  c(9, 3,1,1,6,6,8),
                    status=c(1,NA,1,0,1,1,0),
                    x=     c(0, 2,1,1,1,0,0),
                    z= factor(c('a', 'a', 'a', 'a', 'c', 'c', 'a')))

ftest <- model.frame(fit1, data=test2)
all.equal(levels(ftest$z), levels(test1$z))

# xtest will have one more row than the others, since it does not delete
#  the observation with a missing value for status
xtest <- model.matrix(fit1, data=test2)
dummy <- fit2$x
dummy[,1] <- 0
all.equal(xtest[-2,], dummy, check.attributes=FALSE)
