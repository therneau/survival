options(na.action=na.exclude) # preserve missings
options(contrasts=c('contr.treatment', 'contr.poly')) #ensure constrast type
library(survival)

aeq <- function(x,y, ...) all.equal(as.vector(x), as.vector(y), ...)

# fit1 and fit4 should follow identical iteration paths
fit1 <- survreg(Surv(futime, fustat) ~ age + ecog.ps, ovarian, x=TRUE)
fit4 <- survreg(Surv(log(futime), fustat) ~age + ecog.ps, ovarian,
		dist='extreme')
aeq(fit1$coef, fit4$coef)
aeq(fit1$var, fit4$var)
 
resid(fit1, type='working')
resid(fit1, type='response')
resid(fit1, type='deviance')
resid(fit1, type='dfbeta')
resid(fit1, type='dfbetas')
resid(fit1, type='ldcase')
resid(fit1, type='ldresp')
resid(fit1, type='ldshape')
resid(fit1, type='matrix')

aeq(resid(fit1, type='working'),resid(fit4, type='working'))
#aeq(resid(fit1, type='response'), resid(fit4, type='response'))#should differ
aeq(resid(fit1, type='deviance'), resid(fit4, type='deviance'))
aeq(resid(fit1, type='dfbeta'),   resid(fit4, type='dfbeta'))
aeq(resid(fit1, type='dfbetas'),  resid(fit4, type='dfbetas'))
aeq(resid(fit1, type='ldcase'),   resid(fit4, type='ldcase'))
aeq(resid(fit1, type='ldresp'),   resid(fit4, type='ldresp'))
aeq(resid(fit1, type='ldshape'),  resid(fit4, type='ldshape'))
aeq(resid(fit1, type='matrix'),   resid(fit4, type='matrix'))

# Test suggested by Achim Zieleis: residuals should give a score vector
r1 <-residuals(fit1, type='matrix')
score <- c(as.vector(r1[,c("dg")]) %*% model.matrix(fit1),
           "log(scale)" = sum(r1[,"ds"]))
all(abs(score) < 1e-6)

# repeat this with Gaussian (no transform = different code path)
tfit <-  survreg(Surv(durable, durable>0, type='left') ~age + quant,
                data=tobin, dist='gaussian')
r2 <- residuals(tfit, type='matrix')
score <- c(as.vector(r2[, "dg"]) %*% model.matrix(tfit),
             "log(scale)" = sum(r2[,"ds"]))
all(score < 1e-6)

#
# Some tests of the quantile residuals
#
# These should agree exactly with Ripley and Venables' book
fit1 <- survreg(Surv(time, status) ~ temp, data= imotor)
summary(fit1)

#
# The first prediction has the SE that I think is correct
#  The third is the se found in an early draft of Ripley; fit1 ignoring
#  the variation in scale estimate, except via it's impact on the
#  upper left corner of the inverse information matrix.
# Numbers 1 and 3 differ little for this dataset
#
predict(fit1, data.frame(temp=130), type='uquantile', p=c(.5, .1), se=T)

fit2 <- survreg(Surv(time, status) ~ temp, data=imotor, scale=fit1$scale)
predict(fit2, data.frame(temp=130), type='uquantile', p=c(.5, .1), se=T)

fit3 <- fit2
fit3$var <- fit1$var[1:2,1:2]
predict(fit3, data.frame(temp=130), type='uquantile', p=c(.5, .1), se=T)

pp <- seq(.05, .7, length=40)
xx <- predict(fit1, data.frame(temp=130), type='uquantile', se=T,
	      p=pp)
#matplot(pp, cbind(xx$fit, xx$fit+2*xx$se, xx$fit - 2*xx$se), type='l')


#
# Now try out the various combinations of strata, #predicted, and
#  number of quantiles desired
#
fit1 <- survreg(Surv(time, status) ~ inst + strata(inst) + age + sex, lung)
qq1 <- predict(fit1, type='quantile', p=.3, se=T)
qq2 <- predict(fit1, type='quantile', p=c(.2, .3, .4), se=T)
aeq <- function(x,y) all.equal(as.vector(x), as.vector(y))
aeq(qq1$fit, qq2$fit[,2])
aeq(qq1$se.fit, qq2$se.fit[,2])

qq3 <- predict(fit1, type='quantile', p=c(.2, .3, .4), se=T,
	       newdata= lung[1:5,])
aeq(qq3$fit, qq2$fit[1:5,])

qq4 <- predict(fit1, type='quantile', p=c(.2, .3, .4), se=T, newdata=lung[7,])
aeq(qq4$fit, qq2$fit[7,])

qq5 <- predict(fit1, type='quantile', p=c(.2, .3, .4), se=T, newdata=lung)
aeq(qq2$fit, qq5$fit)
aeq(qq2$se.fit, qq5$se.fit)
