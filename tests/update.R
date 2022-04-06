library(survival)

# the way a +cluster() term is handled in coxph has implications for update.

fit1 <- coxph(Surv(time, status) ~ age, cluster= inst, lung)
fit2 <- coxph(Surv(time, status) ~ age + cluster(inst), lung)
all.equal(fit1, fit2)

fit3 <- coxph(Surv(time, status) ~ age + sex + cluster(inst), lung)

test1 <- update(fit1, .~ .+ sex)
all.equal(test1, fit3)

# Gives a spurious warning message
test2 <- update(fit1, . ~ age + sex + cluster(inst), lung)
all.equal(test2, fit3)

