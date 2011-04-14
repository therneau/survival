options(na.action=na.exclude) # preserve missings
options(contrasts=c('contr.treatment', 'contr.poly')) #ensure constrast type
library(survival)

#
# Test out the survfit.ci function, which does competing risk
#   estimates
#
# For this we need the sequential MGUS data set, using the first
#   obs for each subject
#
aeq <- function(x,y, ...) all.equal(as.vector(x), as.vector(y), ...)

tdata <- data.frame(time=mgus1$stop,
                    status=mgus1$status,
                    event=mgus1$event,
                    sex=mgus1$sex)[mgus1$start==0,]

fit1 <- survfit(Surv(time, status) ~ 1, etype=event, tdata)

# Now get the overall survival, and the hazard for progression
fit2 <- survfit(Surv(time, status) ~1, tdata)  #overall to "first bad thing"
fit3 <- survfit(Surv(time, status*(event=='progression')) ~1, tdata,
                type='fleming')

# Classic CI formula
#  integral [hazard(t) S(t-0) dt], where S= "survival to first event"
haz <- diff(c(0, -log(fit3$surv))) #Aalen hazard estimate
tsurv <- c(1, fit2$surv[-length(fit2$surv)])  #lagged survival
ci <- cumsum(haz *tsurv)
aeq(1-ci, fit1$surv[,1])

#
# Now, make sure that it works for subgroups
#
fit1 <- survfit(Surv(time, status) ~ sex, etype=event, tdata)
fit2 <- survfit(Surv(time, status) ~ 1, etype=event, tdata,
                        subset=(sex=='male'))
fit3 <- survfit(Surv(time, status) ~ 1, etype=event, tdata,
                   subset=(sex=='female'))

aeq(fit2$surv, fit1$surv[1:fit1$strata[1],])
aeq(fit3$surv, fit1$surv[-(1:fit1$strata[1]),])

#  A second test of cumulative incidence
# compare results to Bob Gray's functions
#  The file gray1 is the result of 
#
#    tstat <- ifelse(tdata$status==0, 0, 1+ (tdata$event=='death'))
#    gray1 <- cuminc(tdata$time, tstat)
load("gray1.rda")
plot(gray1[[1]]$time, gray1[[1]]$est, type='l',
     ylim=range(c(gray1[[1]]$est, gray1[[2]]$est)),
     xlab="Time")
lines(gray1[[2]]$time, gray1[[2]]$est, lty=2)

fit2 <- survfit(Surv(time, status) ~ 1, etype=event, tdata)
matlines(fit2$time, 1-fit2$surv, col=2, lty=1:2, type='s')

# To formally match these is a bit of a nuisance.
#  The cuminc function returns full step function, and survfit.ci only
# the bottoms of the steps.
#  The survfit.ci function returns all time points, cuminc only the jumps.
temp1 <- tapply(gray1[[1]]$est, gray1[[1]]$time, max)
indx1 <- match(names(temp1), fit2$time)
aeq(temp1, 1-fit2$surv[indx1,1])
    
