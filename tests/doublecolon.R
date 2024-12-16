# 
# Check that my updates to a. remove survival:: out of formulas and
#  b: ensure that Surv(), cluster(), strata(), pspline(), and tt() use
#   these functions from the survival namespace, when called in a coxph,
#   survfit, survreg, etc formula
library(survival)
aeq <- function(x,y, ...) all.equal(as.vector(x), as.vector(y))

c1 <- coxph(Surv(time, status) ~ age + strata(inst), lung)

# a local Surv that gives the wrong answer but won't error out (makes it
#  simpler to write the tests)
Surv <- function(x, ...) survival::Surv(x, rep(1, length(x)))
c2 <-  coxph(Surv(time, status) ~ age + strata(inst), lung)

c3 <- coxph(survival::Surv(time, status) ~ age + survival::strata(inst), lung)
# in prior releases the above fits a different model, stata is not recognized
#  as a special and becomes a factor

all.equal(coef(c1), coef(c2))
all.equal(coef(c1), coef(c3))

!(c2$call$formula == c3$call$formula)  # c3$call will have 2 survival::, c2 none
deparse1(c3$formula) == "survival::Surv(time, status) ~ age + strata(inst)"  

nocall <- function(x, omit="call") {
    z <- unclass(x)  # needed for any object with a [ method
    z[-match(omit, names(z))]  # all but $call
}
y2 <- with(lung, survival::Surv(time, status)) # outside a formula

fit1a <- coxph(Surv(time, status) ~ age + strata(sex) + cluster(inst), lung)
fit1b <- coxph(Surv(time, status) ~ age + survival::strata(sex) +
                   survival::cluster(inst), lung)
fit1c <- coxph(y2 ~ age + strata(sex) + survival::cluster(inst), lung)
all.equal(nocall(fit1a), nocall(fit1b))
aeq(coef(fit1a), coef(fit1c))

fit2a <- survdiff(Surv(time, status) ~ sex + strata(inst), lung)
fit2b <- survdiff(Surv(time, status) ~ sex + survival::strata(inst),
                  data= lung)
all.equal(nocall(fit2a), nocall(fit2b))
aeq(rowSums(fit2a$obs), c(111, 53))  # make sure it use the correct Surv

fit3a <- survreg(Surv(time, status) ~ ph.ecog + strata(sex), lung)
fit3b <- survreg(Surv(time, status) ~ ph.ecog + survival::strata(sex),
                 data= survival::lung)
all.equal(nocall(fit3a), nocall(fit3b))

fit4a <- concordance(Surv(time, status) ~ ph.ecog + strata(sex), lung)
fit4b <- concordance(Surv(time, status) ~ ph.ecog + survival::strata(sex), lung)
all.equal(nocall(fit4a), nocall(fit4b))

fit5a <- survfit(Surv(time, status) ~ sex, lung)
fit5b <- survfit(Surv(time, status) ~ strata(sex), lung)
fit5c <- survfit(Surv(time, status) ~  survival::strata(sex), lung)
fit5d <- survfit(y2 ~  survival::strata(sex), lung)
all.equal(nocall(fit5a, c("call", "strata")), nocall(fit5b, c("call", "strata")))
all.equal(nocall(fit5b), nocall(fit5c))
aeq(fit5a$surv, fit5d$surv)
