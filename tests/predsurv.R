#
# Verify that predict(coxfit) agrees with survfit for type= expected and survival
#   This also acts as a check on summary.survfit
#
library(survival)
aeq <- function(x,y, ...) all.equal(as.vector(x), as.vector(y), ...)

fit1 <- coxph(Surv(time, status) ~ age + ph.ecog, lung)

d1 <- data.frame(age= 40 * 1:8*5, ph.ecog= rep(0:2, length=8))
curves <- survfit(fit1, newdata= d1)

# the status variable isn't used, but has to be there
d1a <- cbind(time= 1:8 * 60, status=1, d1)  # unique time per subject
d1b <- cbind(time=  365, status=1, d1)      # same time for each

p1 <- predict(fit1, newdata= d1a, type="expected", se.fit=TRUE)
p2 <- predict(fit1, newdata= d1a, type="survival", se.fit=TRUE)
p3 <- predict(fit1, newdata= d1b, type="survival", se.fit=TRUE)

csum1 <- summary(curves, time= 1:8 * 60)
aeq(p1$fit, diag(csum1$cumhaz))
aeq(p1$se.fit , diag(csum1$std.chaz))
aeq(p2$fit, diag(csum1$surv))
aeq(p2$se.fit, diag(csum1$std.err))

csum2 <- summary(curves, time=365)
aeq(p3$fit,csum2$surv)
aeq(p3$se.fit, csum2$std.err)

# Harder, add a strata
fit2 <- coxph(Surv(time, status) ~ age + ph.ecog + strata(sex), lung)
d2 <- data.frame(age= 40 + 1:8*5, ph.ecog= rep(0:2, length=8), sex=rep(1:2,4))
curve2 <- survfit(fit2, newdata= d2[,1:2])

d2a <- cbind(time= 1:8 * 60, status=1, d2)  # unique time per subject
d2b <- cbind(time=  365, status=1, d2)      # same time for each

p1 <- predict(fit2, newdata= d2a, type="expected", se.fit=TRUE)
p2 <- predict(fit2, newdata= d2a, type="survival", se.fit=TRUE)
p3 <- predict(fit2, newdata= d2b, type="survival", se.fit=TRUE)

csum1 <- summary(curve2, time= 1:8 * 60, data.frame=TRUE) 
dummy <- data.frame(d2[,1:2], time= 1:8*60, strata=paste0("sex=", rep(1:2,4)))
temp <- merge(dummy, csum1, all.x=TRUE) # select the correct rows from csum1
aeq(p1$fit, temp$cumhaz) 
aeq(p1$se.fit , temp$std.chaz)
aeq(p2$fit, temp$surv)
aeq(p2$se.fit, temp$std.err)

csum2 <- summary(curve2, time=365)
indx <- cbind(rep(1:2, 4), 1:8)
aeq(p3$fit,csum2$surv[indx])
aeq(p3$se.fit, csum2$std.err[indx])

# Repeat for (time1, time2) data
fit3 <- coxph(Surv(tstart, tstop, status) ~ age + treat, cgd)
d3 <- data.frame(age= c(2,12, 20,30, 40), 
                 treat = rep(levels(cgd$treat), length=5))

curve3 <- survfit(fit3, newdata=d3)

d3a <- cbind(tstart= c(0, 50, 100, 200, 250),
             tstop = c(400, 140, 150, 260, 310), status=1, d3)
d3b <- cbind(tstart=0, tstop=365, status=1, d3)

p4 <- predict(fit3, newdata=d3a, type='expected', se.fit=TRUE)
p5 <- predict(fit3, newdata=d3b, type='survival', se.fit=TRUE)
# type survival is only valid from the start of the curve forward, so no d3a

alltime <- sort(unique(c(d3a$tstart, d3a$tstop)))
csum3 <-summary(curve3, times= alltime)
temp1 <- csum3$cumhaz[cbind(match(d3a$tstart, alltime), 1:nrow(d3a))]
temp2 <- csum3$cumhaz[cbind(match(d3a$tstop, alltime), 1:nrow(d3a))]
aeq(p4$fit, temp2- temp1)

temp3 <- csum3$std.chaz[cbind(match(d3a$tstart, alltime), 1:nrow(d3a))]
temp4 <- csum3$std.chaz[cbind(match(d3a$tstop, alltime), 1:nrow(d3a))]
aeq(p4$se.fit, sqrt(temp4^2 - temp3^2))

csum4 <- summary(curve3, times=365)
aeq(p5$fit, csum4$surv)
aeq(p5$se.fit, csum4$std.err)


# Harder case: add a strata to the problem
fit4 <- coxph(Surv(tstart, tstop, status) ~ age + treat + strata(hos.cat), cgd)
d4 <- data.frame(age= c(2,12, 20,30, 40), 
                 treat = rep(levels(cgd$treat), length=5),
                 hos.cat= levels(cgd$hos.cat)[c(1,2,4,3,2)])

curve4 <- survfit(fit4, newdata=d4)
# by including hos.cat in the above curve4 has 5 strata, one per subject,
#  and no data dimension

d4a <- cbind(tstart= c(0, 50, 100, 200, 250),
             tstop = c(400, 140, 150, 260, 310), status=1, d4)
d4b <- cbind(tstart=0, tstop=365, status=1, d4)
p4 <- predict(fit4, newdata=d4a, type='expected', se.fit=TRUE)
p5 <- predict(fit4, newdata=d4b, type='survival', se.fit=TRUE)

# sort() skipped on purpose, summary should handle it
alltime <- unique(c(d4a$tstart, d4a$tstop))                        
csum5 <- summary(curve4, times=alltime, extend=TRUE)
indx1 <- match(d4a$tstart, csum5$time) + 0:4*length(alltime)
indx2 <- match(d4a$tstop,  csum5$time) + 0:4*length(alltime)
aeq(p4$fit, csum5$cumhaz[indx2] - csum5$cumhaz[indx1])
aeq(p4$se.fit, sqrt(csum5$std.chaz[indx2]^2 - csum5$std.chaz[indx1]^2))

csum6 <- summary(curve4, times= 365, extend=TRUE)
aeq(p5$fit, csum6$surv)
aeq(p5$se.fit, csum6$std.err)
