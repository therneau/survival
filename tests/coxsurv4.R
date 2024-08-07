library(survival)

# Strata by covariate interactions, a case pointed out in early 2011
#  by Frank Harrell, which as it turns out had never been computed
#  correctly by any version of the package.  Which shows how often this
#  case arises in practice.
#
aeq <- function(x, y, ...) all.equal(as.vector(x), as.vector(y))
fit1 <- coxph(Surv(time, status) ~ wt.loss + age*strata(sex) + strata(ph.ecog),
              data=lung)
tdata <- data.frame(wt.loss=c(10,5,0,10, 15,20,25),
                    age    =c(50,60,50,60,70,40,21),
                    sex    =c(1,1,2,2,1,1,1),
                    ph.ecog=c(0,0,1,1,2,2,2))
surv1 <- survfit(fit1, newdata=tdata)

fit2 <- coxph(Surv(time, status) ~ wt.loss + age + I(age*0), data=lung, 
              init=fit1$coefficients, iter=0, subset=(sex==1 & ph.ecog==0))
fit2$var <- fit1$var

surv2 <- survfit(fit2, newdata=list(wt.loss=c(10,5), age=c(50,60)))
s1 <- surv1[1:2]
aeq(s1$surv, surv2$surv) #first a vector, second a matrix
aeq(s1$std.err, surv2$std.err)
aeq(s1[1]$time, surv2$time)
aeq(s1[1]$n.event, surv2$n.event)

fit3 <- coxph(Surv(time, status) ~ wt.loss + age + I(age*1),
              data=lung, init=fit1$coefficients, iter=0,
              subset=(sex==2 & ph.ecog==1))
fit3$var <- fit1$var
surv3 <- survfit(fit3, newdata=list(wt.loss=c(0,10), age=c(50,60)))
aeq(surv1[3:4]$surv, surv3$surv)
aeq(surv1[3:4]$std, surv3$std)

fit4 <- coxph(Surv(time, status) ~ wt.loss + age + I(age*0),
              data=lung, init=fit1$coefficients, iter=0,
              subset=(sex==1 & ph.ecog==2))
fit4$var <- fit1$var
surv4 <- survfit(fit4, newdata=list(wt.loss=c(15,20,25), age=c(70,40,21)))

aeq(surv1[5:7]$surv, surv4$surv)
aeq(surv1[5:7]$std.err, surv4$std.err)
aeq(surv1[5]$n.risk, surv4$n.risk)
             
