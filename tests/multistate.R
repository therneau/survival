#
# Tests for multi-state Cox models
#
library(survival)

aeq <- function(x,y, ...) all.equal(as.vector(x), as.vector(y), ...)

#  There are a few subjects with progression and death on the same day. In the
# usual multi-state data set only one will count. 
data1 <- mgus2
data1$etime <- with(data1, ifelse(pstat==1, ptime, futime))
data1$event <- factor(ifelse(data1$pstat==1, 1, 2L*data1$death),
                      0:2, c("censor", "PCM", "death"))

# direct data set with 2 rows per subject, much like mstate package would do
data2 <- mgus2[rep(1:nrow(mgus2) ,2), c("id", "age", "sex", "mspike")]
data2$time <- rep(data1$etime, 2)
data2$status <- 1* c(data1$event=="PCM", data1$event=="death")
data2$type   <- rep(c(2:3), each=nrow(mgus2))

fit1 <- coxph(Surv(etime, event) ~ age + sex + mspike, data1, id=id, x=TRUE,
              robust=FALSE)
fit1a <- coxph(Surv(etime, event=="PCM") ~ age + sex + mspike, data1)
fit1b <- coxph(Surv(etime, event=='death') ~ age + sex + mspike, data1)
fit1c <- coxph(Surv(time, status) ~ strata(type)/(age + sex+ mspike), 
               data2, x=TRUE)

aeq(fit1$loglik, fit1a$loglik + fit1b$loglik)
aeq(fit1$coef, c(fit1a$coef, fit1b$coef))
aeq(fit1$var[1:3, 1:3], fit1a$var)
aeq(fit1$var[4:6, 4:6], fit1b$var)
aeq(fit1$x[,c(1,4,2,5,3,6)], fit1c$x)
aeq(fit1$coef[c(1,4,2,5,3,6)], fit1c$coef)

# force a common age effect across all states
fit2 <- coxph(list(Surv(etime, event) ~ sex,
                                  1:0 ~ age / common), 
              data1, id=id)

data2 <- rbind(cbind(data1, status= (data1$event=="PCM"), etype=1),
               cbind(data1, status= (data1$event=='death'), etype=2))
fit2a <- coxph(Surv(etime, status) ~ age + strata(etype)/sex, data2)

aeq(coef(fit2), coef(fit2a)[c(2,1,3)])  # not in the same order
aeq(fit2$loglik, fit2a$loglik)

#same fit in more complex ways
data1$entry <- "Entry"
fit2b <-  coxph(list(Surv(etime, event) ~ sex,
                     "Entry":"PCM" + "Entry":"death" ~ age / common),
                istate=entry, data1, id=id)
fit2c <-  coxph(list(Surv(etime, event) ~ sex,
                     "Entry":state(c("PCM", "death")) ~ age / common),
                istate=entry, data1, id=id)

aeq(fit2b$loglik, fit2$loglik)
aeq(fit2c$coef, fit2$coef)

# mspike size as a covariate for PCM only
# first, 4 different ways to write the same
fit3 <- coxph(list(Surv(etime, event) ~ age + sex,
                   1:state("PCM") ~ mspike),
              data1, id=id)
fit3b <- coxph(list(Surv(etime, event) ~ age + sex,
                   1:"PCM" ~ mspike),
              data1, id=id)
fit3c <- coxph(list(Surv(etime, event) ~ age + sex,
                   1:c("PCM") ~ mspike),
              data1, id=id)
fit3d <- coxph(list(Surv(etime, event) ~ age + sex + mspike,
                    1:3 ~ -mspike), data1, id=id)

aeq(fit3b$coef, fit3$coef)
aeq(fit3c$coef, fit3$coef)
aeq(fit3d$coef, fit3$coef)

data3 <- data2
data3$mspike[data3$etype==2] <- 0
fit3a <-  coxph(Surv(etime, status) ~ strata(etype)/(age + sex + mspike), data3)
aeq(fit3$loglik, fit3a$loglik)
aeq(fit3$coef, fit3a$coef[c(1,3,5,2,4)])

# models with strata
test1 <-  coxph(Surv(etime, event=="PCM") ~ age + mspike + strata(sex), data1)
test2 <-  coxph(Surv(etime, event=="death") ~ age + strata(sex), data1)

sfit1 <-  coxph(list(Surv(etime, event) ~ age + strata(sex), 
                   1:state("PCM") ~ mspike),
              data1, id=id)
aeq(coef(sfit1), c(coef(test1), coef(test2)))

test3 <- coxph(Surv(etime, event=="death") ~ age +sex, data1)
sfit2 <- coxph(list(Surv(etime, event) ~ age + sex,
                    1:2 ~ mspike + strata(sex) - sex),  data1, id=id)
aeq(coef(sfit2), c(coef(test1), coef(test3)))

