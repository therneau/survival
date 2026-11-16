library(survival)
aeq <- function(x, y, ...) all.equal(as.vector(x), as.vector(y), ...)
#
# The Peto lower limit uses var(S) = S^2 (1-S)/n.risk.  For the default log
#  transform the lower limit is S * exp(-z * sqrt((1-S)/n.risk)), whether or
#  not the std.err component is a robust (infinitesimal jackknife) one.
#
zz <- qnorm(.975)
fit1 <- survfit(Surv(time, status) ~ 1, data=aml, conf.lower="peto")
fit2 <- survfit(Surv(time, status) ~ 1, data=aml, conf.lower="peto",
                robust=TRUE, id=1:nrow(aml))
c(fit1$logse, fit2$logse)
keep <- (fit1$surv > 0 & fit1$surv < 1)
peto <- with(fit1, surv * exp(-zz * sqrt((1-surv)/n.risk)))
aeq(fit1$lower[keep], peto[keep])
aeq(fit2$lower[keep], peto[keep])

# The upper limit is not modified
aeq(fit2$upper, survfit(Surv(time, status) ~ 1, data=aml, robust=TRUE,
                        id=1:nrow(aml))$upper)

# Non-integer case weights imply a robust variance
wt <- rep(c(.5, 1.5), length=nrow(aml))
fit3 <- survfit(Surv(time, status) ~ x, data=aml, weights=wt, conf.lower="peto")
fit3$logse
keep <- (fit3$surv > 0 & fit3$surv < 1)
peto <- with(fit3, surv * exp(-zz * sqrt((1-surv)/n.risk)))
aeq(fit3$lower[keep], peto[keep])
