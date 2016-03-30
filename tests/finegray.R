library(survival)

# Test data set for Fine-Gray regression
fdata <- data.frame(time  =c(1,2,3,4,4,4,5,5,6,8,8, 9,10,12),
                    status=factor(c(1,2,0,1,0,0,2,1,0,0,2, 0,1 ,0), 0:2,
                             c("cen", "type1", "type2")),      
                    x     =c(5,4,3,1,2,1,1,2,2,4,6,1,2, 0),
                    id = 1:14)
test1 <- finegray(Surv(time, status) ~., fdata)
test2 <- finegray(Surv(time, status) ~x, fdata, event=2)

# When creating the censoring time distribution remember that
#  censors happen after deaths, so the distribution does not drop until
#  time 3+, 4+, 6+, 8+ and 9+
csurv <- list(time=c(0, 3, 4, 6, 8, 9),
              p = cumprod(c(1, 11/12, 8/10, 5/6, 3/4, 2/3)))
#
# For estimation of event type 1, the first subject of event type
#  2 will have weights of curve$p over (0,3], (3,4], (4,6], (6,8], (8,9] 
#  and (9,12].  All that really matters is the weight at times 1, 4, 5,
#  and 10, however, which are the points at which events of type 1 happen
#  
# The next subject of event type 2 occurs at time 5, and will have a
#  weight of (9,12] /(4,5] = (5*4*2)/(7*5*3) = 8/21 at time 10.  The last
#  censor at time 6 has a weight of 2/3 at time 10.

tid <- c(1, 2,2,2,2, 3:5, 6,6, 7:10, 11, 11, 12:14)
temp1 <- data.frame(id=tid, start=0, fdata[tid,],
                    wt= c(1, csurv$p[c(1,2,3,6)], 1,1,1,1, 8/21, 1,1,1, 1, 
                         1, 2/3, 1,1,1))
temp1$time <- [c(2,3,4,5, 9,10, 15, 16)] <- c(3,4,5,10, 5,10, 5,10)
temp1$start[c(2,3,4,5, 9,10, 15,16)] <- c(0,3,4,5, 0,5, 0,5)
row.names(temp1) <- NULL
tfit <- survfit(Surv(start, time, status==1) ~1, temp1, weight=wt)


fdata$stat2 <- as.numeric(fdata$status) -1
ctest <- crprep("time", "stat2", data=fdata, trans=1:2, cens=0, 
                keep=c("id", 'x'))
subset(ctest, failcode==1)

xfit <- survfit(Surv(Tstart, Tstop, status==1) ~1, ctest, weight=weight.cens,
                subset=(failcode==1))

