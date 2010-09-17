options(na.action=na.exclude) # preserve missings
options(contrasts=c('contr.treatment', 'contr.poly')) #ensure constrast type
library(survival)

#
# Verify that scale can be fixed at a value
#    coefs will differ slightly due to different iteration paths
tol <- .001

# Intercept only models
fit1 <- survreg(Surv(time,status) ~ 1, lung)
fit2 <- survreg(Surv(time,status) ~ 1, lung, scale=fit1$scale)
all.equal(fit1$coef, fit2$coef, tolerance= tol)
all.equal(fit1$loglik, fit2$loglik, tolerance= tol)

# multiple covariates
fit1 <- survreg(Surv(time,status) ~ age + ph.karno, lung)
fit2 <- survreg(Surv(time,status) ~ age + ph.karno, lung,
		scale=fit1$scale)
all.equal(fit1$coef, fit2$coef, tolerance=tol)
all.equal(fit1$loglik[2], fit2$loglik[2], tolerance=tol)

# penalized models
fit1 <- survreg(Surv(time, status) ~ pspline(age), lung)
fit2 <- survreg(Surv(time, status) ~ pspline(age), lung, scale=fit1$scale)
all.equal(fit1$coef, fit2$coef, tolerance=tol)
all.equal(fit1$loglik[2], fit2$loglik[2], tolerance=tol)

rm(fit1, fit2, tol)

