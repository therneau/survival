library(survival)

# Strata by covariate interactions, a case pointed out in early 2011
#  by Frank Harrell, which as it turns out had never been computed
#  correctly by any version of the package.  This shows how often this
#  case arises in practice.
#

fit1 <- coxph(Surv(time, status) ~ wt.loss + age*strata(sex) + strata(ph.ecog),
              data=lung)
surv1 <- survfit(fit1, newdata=list(wt.loss=c(10,5), age=c(50,60)))

fit2 <- coxph(Surv(time, status) ~ wt.loss + age + I(age*0), data=lung, 
              init=fit1$coef, iter=0, subset=(sex==1 & ph.ecog==0))
fit2$var <- fit1$var

surv2 <- survfit(fit2, newdata=list(wt.loss=c(10,5), age=c(50,60)))

nocall <- length(surv2)  # the calls won't be the same
all.equal(unclass(surv1[1,])[-nocall], unclass(surv2)[-nocall])

fit3 <- coxph(Surv(time, status) ~ wt.loss + age + I(age*0),
              data=lung, init=fit1$coef, iter=0,
              subset=(sex==1 & ph.ecog==1))
fit3$var <- fit1$var
surv3 <- survfit(fit3, newdata=list(wt.loss=c(10,5), age=c(50,60)))
all.equal(surv1[2,]$surv, surv3$surv)
all.equal(surv1[2,]$std, surv3$std)

fit4 <- coxph(Surv(time, status) ~ wt.loss + age + I(age),
              data=lung, init=fit1$coef, iter=0,
              subset=(sex==2 & ph.ecog==1))
fit4$var <- fit1$var
surv4 <- survfit(fit4, newdata=list(wt.loss=c(10,5), age=c(50,60)))

all.equal(surv1[6,]$surv, surv4$surv)
all.equal(surv1[6,]$std.err, surv4$std.err)
             
fit5 <- coxph(Surv(time, status) ~ wt.loss + age + I(age) + strata(ph.ecog),
              data=lung, init=fit1$coef, iter=0,
              subset=(sex==2))
fit5$var <- fit1$var
surv5 <- survfit(fit5, newdata=list(wt.loss=c(10,5), age=c(50,60)))
temp1 <- surv1[5:7,]
temp2 <- surv5[1:3,]
all.equal(temp1$strata, temp2$strata, check.attributes=FALSE)
nlist <- c('n', 'time', 'n.risk', 'n.event', 'n.censor', 'surv', 
           'std.err')
all.equal(unclass(temp1)[nlist], unclass(temp2)[nlist])
