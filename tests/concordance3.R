library(survival)
options(na.action=na.exclude)
aeq <- function(x,y,...)  all.equal(as.vector(x), as.vector(y), ...)

#  Make sure strata is retained, and that the overall variance is correct
fit1 <- coxph(Surv(time, status) ~ age + offset(ph.ecog*0) +strata(sex), lung)
fit2 <- coxph(Surv(time, status) ~ age + ph.ecog +strata(sex), lung)

test <- concordance(fit1, fit2, influence=1)

ksex <- model.frame(fit1)[["strata(sex)"]]
test1 <- concordance(fit1$y ~ fit1$linear.predictors + strata(ksex), 
                     reverse=TRUE, influence=1)
test2 <- concordance(fit1$y ~ fit2$linear.predictors + strata(ksex), 
                     reverse=TRUE, influence=1)
aeq(test$concordance, c(test1$concordance, test2$concordance))
aeq(diag(test$var), c(test1$var[1], test2$var[1]))
aeq(test$dfbeta, cbind(test1$dfbeta, test2$dfbeta))

cvec <- c(-1, 1)
aeq(cvec %*% test$var %*% cvec, sum((test1$dfbeta - test2$dfbeta)^2))

# Time weights
# The mgus2 data set has very long follow-up and non-proportional hazards
#
mfit <- coxph(Surv(futime, death) ~ creat + hgb, mgus2)
cm1 <- concordance(mfit, timewt='n', ranks=TRUE)
cm2 <- concordance(mfit, timewt='S', ranks=TRUE)
cm3 <- concordance(mfit, timewt='S/G', ranks=TRUE)
sfit <- survfit(Surv(futime, death) ~ 1, mgus2, subset=!is.na(creat+hgb))
gfit <- survfit(Surv(futime, 1-death)~1, mgus2, subset=!is.na(creat+hgb))

rd1 <- cm1$ranks
rd2 <- cm2$ranks
rd3 <- cm3$ranks
all.equal(rd1[c('time', 'rank', 'casewt')], rd2[c('time', 'rank', 'casewt')])
all.equal(rd1[c('time', 'rank', 'casewt')], rd3[c('time', 'rank', 'casewt')])

# n(t) is the number of comparable pairs, which does not count tied times
indx <- match(rd1$time, sfit$time)
nt <- sfit$n.risk[indx] - sfit$n.event[indx]
all.equal(rd1$timewt, nt)

# I need S(t) and G(t-), G has the exact same time points as S though
gminus <- c(1, gfit$surv)
all.equal(rd2$timewt, mfit$n* sfit$surv[indx])
all.equal(rd3$timewt, mfit$n* sfit$surv[indx] /gminus[indx])


