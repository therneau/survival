#
# Verify that using multiple states + proportional baselines
#  will mimic a factor covariate
#
library(survival)
aeq <- function(x, y, ...) all.equal(as.vector(x), as.vector(y), ...)

tdata <- subset(lung, ph.ecog < 3)  # there is only one of them
tdata$state <- factor(tdata$status, 1:2, c("censor", "death"))
tdata$cstate<- factor(tdata$ph.ecog, 0:2, c("ph0", "ph1", "ph2"))
tdata$id  <- 1:nrow(tdata)
survcheck(Surv(time, state) ~ 1, id=id, istate=cstate, tdata)

fit1 <- coxph(Surv(time, status) ~ age + sex + factor(ph.ecog), tdata, 
                 ties="breslow")
fit2 <- coxph(list(Surv(time, state) ~1,
                   1:4 + 2:4 + 3:4~ age + sex/ common + shared), 
              id=id, istate=cstate, data= tdata, ties="breslow")
aeq(coef(fit1), coef(fit2))
all.equal(fit1$loglik, fit2$loglik)

csurv1 <- survfit(fit1, newdata=expand.grid(age=65, sex=1, ph.ecog=0:2))
csurv2a <- survfit(fit2, newdata= list(age=65, sex=1), p0=c(1,0,0,0))
csurv2b <- survfit(fit2, newdata= list(age=65, sex=1), p0=c(0,1,0,0))
csurv2c <- survfit(fit2, newdata= list(age=65, sex=1), p0=c(0,0,1,0))

aeq(csurv1[1]$surv, csurv2a$pstate[,1,1])
aeq(csurv1[2]$surv, csurv2b$pstate[,1,2])
aeq(csurv1[3]$surv, csurv2c$pstate[,1,3])

# Since the multi-state does not yet implement the Efron approx (and may never)
#  for the survival curve, the above only perfectly matches with Breslow.
