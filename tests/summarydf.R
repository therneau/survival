#
# Check out the data.frame option of summary.survfit
#
library(survival)
aeq <- function(x, y, ...) all.equal(as.vector(x), as.vector(y), ...)

surv1 <- survfit(Surv(time, status) ~ ph.ecog, lung)
summ1a <- summary(surv1, times= round(30.5*1:32))
summ1b <- summary(surv1,  times= round(30.5*1:32), data.frame=TRUE)
test <- with(summ1a, data.frame(time, n.risk, n.event, n.censor, 
              surv, cumhaz, std.err, std.chaz, lower, upper, strata))
all.equal(test, summ1b)

# more complex
cox1 <- coxph(Surv(time, status) ~ ph.ecog + sex + strata(inst), lung)
dummy <- expand.grid(ph.ecog=0:2, sex=1:2)
surv2 <- survfit(cox1, newdata=dummy)

summ2a <- summary(surv2, times= round(30.5*1:32)) 
summ2b <- summary(surv2, times= round(30.5*1:32), data.frame=TRUE)
nd <- 6 # the "data" dimension
ns <- length(summ2a$n.risk)

# The result should have all 19 strata for newdata[1,], then for newdata[2,]
test2 <- with(summ2a, data.frame(time=rep(time, nd), n.risk=rep(n.risk, nd),
                                 n.event = rep(n.event, nd),
                                 n.censor= rep(n.censor, nd),
                                 surv= c(surv), cumhaz=c(cumhaz),
                                 std.err=c(std.err), std.chaz=c(std.chaz),
                                 lower=c(lower), upper=c(upper), 
                                 strata=rep(strata, nd),
                                 dummy[rep(1:nd, each=ns),]))
all.equal(test2, summ2b)

# Now a multistate one
mgus2$etime <- with(mgus2, ifelse(pstat==0, futime, ptime))
temp <- with(mgus2, ifelse(pstat==0, 2*death, 1))
mgus2$event <- factor(temp, 0:2, labels=c("censor", "pcm", "death"))

surv3 <- survfit(Surv(etime,event) ~ sex, mgus2, id=id)
summ3a<- summary(surv3, times= 1:30 *12)  # yearly
summ3b <-summary(surv3, times =1:30 *12, data.frame=TRUE)

# 3 states, results for state 1, then state 2, etc.
test3 <- with(summ3a, data.frame(time= rep(time,3), n.risk = c(n.risk),
                                 n.event = c(n.event),
                                 n.censor= c(n.censor), pstate=c(pstate),
                                 std.err= c(std.err), lower=c(lower),
                                 upper=c(upper), strata=rep(strata,3),
                           state= rep(surv3$state, each=length(summ3a$time))))
all.equal(test3, summ3b)

cox4 <- coxph(Surv(etime, event) ~ age + mspike + strata(sex), mgus2, id=id)
d4 <- expand.grid(age=c(70, 80, 90), mspike=c(0, 1))

surv4 <- survfit(cox4, newdata=d4)
dim4 <- dim(surv4)
dim4

# The data frame will be consistent with R matrices, i.e., the first
#  subscript varies fastest: strata, then data, then states
summ4a <- summary(surv4, times= 1:30 *12)
summ4b <- summary(surv4, times= 1:30 *12, data.frame=TRUE)

j <- rep(1:length(summ4a$time), dim4[2])
k <- rep(rep(1:dim4[2], each= length(summ4a$time)), dim4[3])
test4 <- with(summ4a, data.frame(time= rep(time,18), n.risk = c(n.risk[j,]),
                                 n.event = c(n.event[j,]),
                                 n.censor= c(n.censor[j,]), pstate=c(pstate),
                                 strata=rep(strata, 18),
                                 state= rep(surv4$state, each=6* length(time)),
                                 d4[k,]))
all.equal(summ4b, test4)

# without strata
cox5 <-coxph(Surv(etime, event) ~ age + mspike, mgus2, id=id) 
surv5 <- survfit(cox5, newdata=d4)
summ5a <- summary(surv5, times= 1:30 *12)
summ5b <- summary(surv5, times= 1:30 *12, data.frame=TRUE)

j <- rep(1:length(summ5a$time), dim4[2])
k <- rep(rep(1:dim4[2], each= length(summ5a$time)), dim4[3])
test5 <- with(summ5a, data.frame(time= rep(time,18), n.risk = c(n.risk[j,]),
                                 n.event = c(n.event[j,]),
                                 n.censor= c(n.censor[j,]), pstate=c(pstate),
                                 state= rep(surv5$state, each=6* length(time)),
                                 d4[k,]))
all.equal(summ5b, test5)

# without times
summ6a <- summary(surv5)
summ6b <- summary(surv5, data.frame=TRUE)
j <- rep(1:length(summ6a$time), dim4[2])
k <- rep(rep(1:dim4[2], each= length(summ6a$time)), dim4[3])
test6 <- with(summ6a, data.frame(time= rep(time,18), n.risk = c(n.risk[j,]),
                                 n.event = c(n.event[j,]),
                                 n.censor= c(n.censor[j,]), pstate=c(pstate),
                                 state= rep(surv5$state, each=6* length(time)),
                                 d4[k,]))
all.equal(summ6b, test6)
