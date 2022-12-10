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
# Start with a very small data set: aml has 23 subjects
#
atest1 <- concordance(Surv(time, status) ~ x, aml, ranks=TRUE)
atest2 <- concordance(Surv(time, status) ~ x, aml, ranks=TRUE, timewt='S')
atest3 <- concordance(Surv(time, status) ~ x, aml, ranks=TRUE, timewt='S/G')
atest4 <- concordance(Surv(time, status) ~ x, aml, ranks=TRUE, timewt='n/G2')
# The ranks data frame agrees for all but weights
all.equal(atest1$ranks[, -3], atest2$ranks[, -3])
all.equal(atest1$ranks[, -3], atest3$ranks[, -3])
all.equal(atest1$ranks[, -3], atest4$ranks[, -3])

wt1 <- cbind(atest1$ranks[,"timewt"], atest2$ranks[,"timewt"],
             atest3$ranks[,"timewt"], atest4$ranks[,"timewt"])
            
# survfit0 adds time 0 to the curves
# to break ties between censor/death for G, we need to add an offset to
#  the censoring times.  Since time is integer, .1 works nicely
s1 <- survfit0(survfit(Surv(time, status) ~ 1, aml))
g1 <- survfit0(survfit(Surv(time + .1*(1-status), 1-status) ~1, aml))

# The ingredients of the weights
indx <- match(atest1$ranks[,"time"], s1$time)
nrisk  <- s1$n.risk[indx]
sminus <- s1$surv[indx-1]
gminus <- g1$surv[findInterval(atest1$ranks[,"time"], g1$time)]
n <- nrow(aml)

wt2 <- cbind(nrisk, n*sminus, n*sminus/gminus, nrisk/gminus^2)
aeq(wt1, wt2)

# The sum of weighted ranks should equal (C-D) for a Cox model fit
tfun <- function(cfit, reverse=FALSE) {
    t1 <- sum(cfit$ranks$timewt * cfit$ranks$rank)
    t2 <- cfit$count[1] - cfit$count[2]
    all.equal(unname(t1), unname(t2))
}
tfun(atest1)
tfun(atest2)
tfun(atest3)
tfun(atest4)

# The nafld data set has strong and early censoring (one of the only ones
# in the package that does.) So it is a good check of time weights.
#
nfit <- coxph(Surv(futime, status) ~ male + pspline(age), nafld1)
cn1 <- concordance(nfit, timewt='n', ranks=TRUE)
cn2 <- concordance(nfit, timewt='S', ranks=TRUE)
cn3 <- concordance(nfit, timewt='S/G', ranks=TRUE)
cn4 <- concordance(nfit, timewt='n/G2', ranks=TRUE)

sfit <- survfit0(survfit(Surv(futime, status) ~ 1, nafld1))
gfit <- survfit0(survfit(Surv(futime + .1*(status==0), 1-status) ~0, nafld1))

# The ingredients of the weights
dtime <- cn1$ranks[, "time"]
indx <- match(dtime, sfit$time)
nrisk  <- sfit$n.risk[indx]
sminus <- sfit$surv[indx-1]
gminus <- gfit$surv[findInterval(dtime, gfit$time)]
n <- nrow(nafld1)

wt1 <- cbind(cn1$ranks[, "timewt"], cn2$ranks[,"timewt"],
             cn3$ranks[, "timewt"], cn4$ranks[,"timewt"])
wt2 <- cbind(nrisk, n*sminus, n*sminus/gminus, nrisk/gminus^2)
aeq(wt1, wt2)

rd1 <- cn1$ranks
rd2 <- cn2$ranks
rd3 <- cn3$ranks
all.equal(rd1[c('time', 'rank', 'casewt')], rd2[c('time', 'rank', 'casewt')])
all.equal(rd1[c('time', 'rank', 'casewt')], rd3[c('time', 'rank', 'casewt')])

tfun(cn1)
tfun(cn2)
tfun(cn3)
tfun(cn4)

