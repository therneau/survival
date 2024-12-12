# 
# Check that my updates to a. remove survival:: out of formulas and
#  b: ensure that Surv(), cluster(), strata(), pspline(), and tt() use
#   those functions from the survival namespace, when called in a coxph,
#   survfit, survreg, etc formula
library(survival)

c1 <- coxph(Surv(time, status) ~ age + strata(inst), lung)

# a local Surv that gives the wrong answer but won't error out
Surv <- function(x, ...) survival::Surv(x, rep(1, length(x)))
c2 <-  coxph(Surv(time, status) ~ age + strata(inst), lung)

c3 <- coxph(survival::Surv(time, status) ~ age + survival::strata(inst), lung)
# in prior releases the above fits a different model, stata is not recognized
#  as a special and becomes a factor

all.equal(coef(c1), coef(c2))
all.equal(coef(c1), coef(c3))

!(c2$call$formula == c3$call$formula)  # the call doesn't get changed
c2$formula == c3$formula               # but the formula is changed


nocall <- function(x) x[-match("call", names(x))]  # all but $call
fit1a <- coxph(Surv(time, status) ~ age + strata(sex) + cluster(inst), lung)
fit1b <- coxph(Surv(time, status) ~ age + survival::strata(sex) +
                   survival::cluster(inst), lung)
all.equal(nocall(fit1a), nocall(fit1b))

fit2a <- survdiff(Surv(time, status) ~ ph.ecog + strata(sex), lung)
fit2b <- survdiff(Surv(time, status) ~ ph.ecog + survival::strata(sex), lung)
all.equal(nocall(fit2a), nocall(fit2b))

fit3a <- survreg(Surv(time, status) ~ ph.ecog + strata(sex), lung)
fit3b <- survreg(Surv(time, status) ~ ph.ecog + survival::strata(sex), lung)
all.equal(nocall(fit3a), nocall(fit3b))

fit4a <- concordance(Surv(time, status) ~ ph.ecog + strata(sex), lung)
fit4b <- concordance(Surv(time, status) ~ ph.ecog + survival::strata(sex), lung)
all.equal(nocall(fit4a), nocall(fit4b))
