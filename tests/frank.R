library(survival)
#
# Check out intercept/interaction for Frank H
#
age2 <- lung$age - 50
fit1 <- coxph(Surv(time, status) ~ age * strata(sex), lung)
fit2 <- coxph(Surv(time, status) ~ age2*strata(sex), lung)

surv1 <- survfit(fit1)
surv2 <- survfit(fit2)
# The call won't match, but the rest should
icall <- match("call", names(surv1))
all.equal(unclass(surv1)[-icall], unclass(surv2)[-icall])


# It should match what I get with a single strata fit

fit3 <- coxph(Surv(time, status) ~ age, data=lung,
              init=fit1$coef[1], subset=(sex==1), iter=0)
surv1b <- survfit(fit3, newdata=fit1$means[1])
icall <- match("call", surv1b)
all.equal(unlist(surv1[1])[-icall], unlist(surv1b)[-icall])

