# check for acceptance of non-syntatic variable names
# survfit was the one with issues
library(survival)
aeq <- function(x,y,...) all.equal(as.vector(x), as.vector(y))

lung2 <- lung
lung2["bad sex"] <- lung2$sex

fit1a <- coxph(Surv(time, status) ~ ph.ecog + sex, lung2)
fit1b <- coxph(Surv(time, status) ~ ph.ecog + `bad sex`, lung2)
aeq(coef(fit1a), coef(fit1b))
aeq(fit1a$loglik, fit1b$loglik)

dummy <- expand.grid(ph.ecog=0:2, sex=1:2)
dummy["bad sex"] <- dummy$sex
surv1a <- survfit(fit1a, newdata= dummy)
surv2a <- survfit(fit1b, newdata= dummy)
all.equal(surv1a$surv, surv2a$surv)


fit2a <- survfit(Surv(time, status) ~ ph.ecog + sex, lung2)
fit2b <- survfit(Surv(time, status) ~ ph.ecog + `bad sex`, lung2)
aeq(fit2a$surv, fit2b$surv)

fit3a <- survreg(Surv(time, status) ~ ph.ecog + sex, lung2)
fit3b <- survreg(Surv(time, status) ~ ph.ecog + `bad sex`, lung2)
aeq(coef(fit3a), coef(fit3b))
aeq(fit3a$loglik, fit3b$loglik)

# Need to add other calls?
