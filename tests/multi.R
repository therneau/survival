#
# Tests for multi-state Cox models
#
library(survival)
aeq <- function(x,y, ...) all.equal(as.vector(x), as.vector(y), ...)

data1 <- mgus2
data1$etime <- with(data1, ifelse(pstat==1, ptime, futime))
data1$event <- factor(ifelse(data1$pstat==1, 1, 2L*data1$death),
                      0:2, c("censor", "PCM", "death"))

fit1 <- coxph(Surv(etime, event) ~ age + sex + mspike, data1, id=id)

fit1a <- coxph(Surv(etime, event=="PCM") ~ age + sex + mspike, data1)
fit1b <- coxph(Surv(etime, event=='death') ~ age + sex + mspike, data1)

aeq(fit1$loglik, fit1a$loglik + fit1b$loglik)
aeq(fit1$coef, c(fit1a$coef, fit1b$coef))
aeq(fit1$var[1:3, 1:3], fit1a$var)
aeq(fit1$var[4:6, 4:6], fit1b$var)

# force a common age effect across all states
fit2 <- coxph(Surv(etime, event) ~ age + sex, data1, id=id,
              covariates=list( 1:1 ~ age / common))

data2 <- rbind(cbind(data1, status= (data1$event=="PCM"), etype=1),
               cbind(data1, status= (data1$event=='death'), etype=2))
fit2a <- coxph(Surv(etime, status) ~ age + strata(etype)/sex, data2)

aeq(coef(fit2), coef(fit2a)[c(2,3,1)])  # not in the same order
aeq(fit2$loglik, fit2a$loglik)

# mspike size as a covariate for PCM only
# first, 3 different ways to write the same
fit3 <- coxph(Surv(etime, event) ~ age + sex, data1, id=id,
              covariates=list(1:state("PCM") ~ mspike))
fit3b <- coxph(Surv(etime, event) ~ age + sex, data1, id=id,
              covariates=list(1:"PCM" ~ mspike))
fit3c <- coxph(Surv(etime, event) ~ age + sex, data1, id=id,
              covariates=list(1:c("PCM") ~ mspike))
aeq(fit3b$coef, fit3$coef)
aeq(fit3c$coef, fit3$coef)

data3 <- data2
data3$mspike[data3$etype==2] <- 0
fit3a <-  coxph(Surv(etime, status) ~ strata(etype)/(age + sex + mspike), data3)
aeq(fit3$loglik, fit3a$loglik)
aeq(fit3$coef, fit3a$coef[c(1,3,2,4,5)])
