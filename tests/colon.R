#
# Verify the success of my "stop the :: madness" fixes
#
library(survival)
aeq <- function(x, y, ...) all.equal(as.vector(x), as.vector(y), ...)

fit1a <- coxph(Surv(time, status) ~ age + strata(sex) + cluster(inst), lung)
fit1b <- coxph(Surv(time, status) ~ age + survival::strata(sex) +
                   survival::cluster(inst), lung)
all.equal(fit1a, fit1b)

fit2a <- survdiff(Surv(time, status) ~ ph.ecog + strata(sex), lung)
fit2b <- survdiff(Surv(time, status) ~ ph.ecog + survival::strata(sex), lung)
all.equal(fit2a, fit2b)

fit3a <- survreg(Surv(time, status) ~ ph.ecog + strata(sex), lung)
fit3b <- survreg(Surv(time, status) ~ ph.ecog + survival::strata(sex), lung)
all.equal(fit3a, fit3b)

fit4a <- concordance(Surv(time, status) ~ ph.ecog + strata(sex), lung)
fit4b <- concordance(Surv(time, status) ~ ph.ecog + survival::strata(sex), lung)
all.equal(fit4a, fit4b)
