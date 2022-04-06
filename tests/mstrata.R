#
# Verify that using multiple states + proportional baselines
#  will mimic a factor covariate
#
library(survival)
aeq <- function(x, y, ...) all.equal(as.vector(x), as.vector(y), ...)

tdata <- subset(lung, ph.ecog < 3)  # there is only one row with ph.ecog=3
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


# part 2: predicted values for a multistate model that has a strata
mgus2$etime <- with(mgus2, ifelse(pstat==0, futime, ptime))
event <- with(mgus2, ifelse(pstat==0, 2*death, 1))
mgus2$event <- factor(event, 0:2, labels=c("censor", "pcm", "death"))

dummy <- expand.grid(sex=c("F", "M"), age=c(60, 80), mspike=1.2)

cfit1 <- coxph(Surv(etime, event) ~ age + mspike +strata(sex), mgus2, id=id)

csurv1 <- survfit(cfit1, newdata=dummy)
csurv2 <- survfit(cfit1, newdata=expand.grid(age=c(60, 80), mspike=1.2))
temp <- c('n', 'time', 'n.risk', 'n.event', 'n.censor', 'pstate', 'cumhaz')
all.equal(unclass(csurv1)[temp], unclass(csurv2)[temp])

cfit2 <- coxph(Surv(etime, event) ~ age + mspike, id=id,
               init= coef(cfit1), iter=0, data=mgus2, subset=(sex=='F'))
csurv3 <- survfit(cfit2, newdata= expand.grid(age=c(60, 80), mspike=1.2))
all.equal(unclass(csurv1[1,,])[temp], unclass(csurv3)[temp])
