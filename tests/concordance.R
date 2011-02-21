options(na.action=na.exclude) # preserve missings
options(contrasts=c('contr.treatment', 'contr.poly')) #ensure constrast type
library(survival)

#
# Simple tests of concordance.  These numbers were derived in multiple
#   codes.
#
aeq <- function(x,y, ...) all.equal(as.vector(x), as.vector(y), ...)

tdata <- aml[aml$x=='Maintained',]
tdata$y <- c(1,6,2,7,3,7,3,8,4,4,5)
fit <- survConcordance(Surv(time, status) ~y, tdata)
aeq(fit$stats, c(14,24,2,0))

# Lots of ties
tempx <- Surv(c(1,2,2,2,3,4,4,4,5,2), c(1,0,1,0,1,0,1,1,0,1))
tempy <- c(5,5,4,4,3,3,7,6,5,4)
fit2 <- survConcordance(tempx ~ tempy)
aeq(fit2$stats, c(13,13,5,2))

# Bigger data
fit3 <- survConcordance(Surv(time, status) ~ age, lung)
aeq(fit3$stats, c(10717, 8706, 591, 28))

# More ties
fit4 <- survConcordance(Surv(time, status) ~ ph.ecog, lung)
aeq(fit4$stats, c(8392, 4258, 7137, 28))

# Case weights
wt <- c(1,2,3,2,1,2,3,4,3,2,1)
fit5 <- survConcordance(Surv(time, status) ~ y, tdata, weight=wt)
fit6 <- survConcordance(Surv(time, status) ~y, tdata[rep(1:11,wt),])
aeq(fit5$stats, c(70, 91, 7, 0))  # checked by hand
aeq(fit5$stats[1:3], fit6$stats[1:3])  #spurious "tied on time" is produced

# Start, stop simplest case
fit7 <- survConcordance(Surv(rep(0,11), time, status) ~ y, tdata, weight=wt)
aeq(fit5$stasts, fit7$stats)

# Multiple intervals for some, but same risk sets as tdata
tdata2 <- data.frame(time1=c(0,3, 5,  6,7,   0,  4,18,  7,  0,27,  2,  0, 
                             0,9, 5),
                     time2=c(3,9, 13, 7,13, 18, 18,23, 28, 27,31, 34, 45, 
                             9,48, 60),
                     status=c(0,1, 1, 0,0,  1,  0,1, 0, 0,1, 1, 0, 0,1, 1),
                     y = c(1,1, 6, 2,2, 7, 3,3, 7, 3,3, 8, 4, 4,4, 5),
                     wt= c(1,1, 2, 3,3, 2, 1,1, 2, 3,3, 4, 3, 2,2, 1))
fit8 <- survConcordance(Surv(time1, time2, status) ~y, tdata2, weight=wt)
aeq(fit5$stats, fit8$stats)


# Stratified
tdata3 <- data.frame(time1=c(tdata2$time1, rep(0, nrow(lung))),
                     time2=c(tdata2$time2, lung$time),
                     status = c(tdata2$status, lung$status -1),
                     x = c(tdata2$y, lung$ph.ecog),
                     wt= c(tdata2$wt, rep(1, nrow(lung))),
                     grp=rep(1:2, c(nrow(tdata2), nrow(lung))))
fit9 <- survConcordance(Surv(time1, time2, status) ~x + strata(grp),
                        data=tdata3, weight=wt)
aeq(fit9$stats[1,], fit5$stats)
aeq(fit9$stats[2,], fit4$stats)

                             
