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
# The mgus2 data set has very long follow-up and non-proportional hazards;
#  so one could argue that time weights are appropriate in this case
#
mfit <- coxph(Surv(futime, death) ~ creat + hgb, mgus2)
cm1 <- concordance(mfit, timewt='n', ranks=TRUE)
cm2 <- concordance(mfit, timewt='S', ranks=TRUE)
cm3 <- concordance(mfit, timewt='S/G', ranks=TRUE)
sfit <- survfit(Surv(futime, death) ~ 1, mgus2, subset=!is.na(creat+hgb))
sfit <- survfit0(sfit)  # at time 0
# To be completely correct, you can't just reverse the status in order to
#  get the censoring function G(t).  Censorings happen after deaths, so ties
#  work differently.
# In the mgus2 data set, there is 1 death and 10 censors at time 10, 1202 at
#  risk at time 10: the KM multiplier for S(10) is 1201/1202.  But the
#  multiplier for G(10) is 1191/1201 -- the death is not "at risk" for 
#  censoring at time 10.
# The futime variable in MGUS is an interger, so add .1 to all the censoring
#  times to move them into the future.
gtime <- ifelse(mgus2$death==1, mgus2$futime, mgus2$futime + .1)
gfit <- survfit(Surv(gtime, 1-death)~1, mgus2, subset=!is.na(creat+hgb))
gfit <- survfit0(gfit)  # add time 0

rd1 <- cm1$ranks
rd2 <- cm2$ranks
rd3 <- cm3$ranks
all.equal(rd1[c('time', 'rank', 'casewt')], rd2[c('time', 'rank', 'casewt')])
all.equal(rd1[c('time', 'rank', 'casewt')], rd3[c('time', 'rank', 'casewt')])

# n(t) is the number of comparable pairs, which does not count tied times
sindx <- findInterval(rd1$time, sfit$time, left.open=TRUE)
gindx <- findInterval(rd1$time, gfit$time, left.open=TRUE) 
all.equal(rd1$timewt, sfit$n.risk[sindx])

nt <- sfit$n.risk[sindx] 
all.equal(rd1$timewt, nt)

# I need S(t-) and G(t-),  at the death unique death times in S
sminus <- c(1, sfit$surv)[1L +findInterval(rd1$time, sfit$time, left.open=TRUE)]
gminus <- c(1, gfit$surv)[1L +findInterval(rd1$time, gfit$time, left.open=TRUE)]
all.equal(rd2$timewt, mfit$n* sminus)
all.equal(rd3$timewt, mfit$n* sminus/gminus)


