options(na.action=na.exclude) # preserve missings
options(contrasts=c('contr.treatment', 'contr.poly')) #ensure constrast type
library(survival)

#
# Test out the t-distribution
#
# First, a t-dist with 500 df should be nearly identical to the Gaussian

fitig <- survreg(Surv(time, status)~voltage, 
        dist = "gaussian", data = capacitor)
fit1 <- survreg(Surv(time, status) ~ voltage,
		 dist='t', parms=500, capacitor)
fitig
summary(fit1, corr=F)

# A more realistic fit
fit2 <- survreg(Surv(time, status) ~ voltage,
		 dist='t', parms=5, capacitor)
print(fit2)

if (FALSE) {
resid(fit2, type='response')
resid(fit2, type='deviance')
resid(fit2, type='working') 
resid(fit2, type='dfbeta')
resid(fit2, type='dfbetas')
resid(fit2, type='ldresp')
resid(fit2, type='ldshape')
resid(fit2, type='ldcase')
resid(fit2, type='matrix')

predict(fit2, type='link')
predict(fit2, type='terms')
predict(fit2, type='quantile')
}
