library(survival)

# Test data set for Fine-Gray regression
fdata <- data.frame(time  =c(1,2,3,4,4,5,5,5,6,7,8,9,10),
                    status=c(1,2,0,1,0,0,0,1,0,2,2,0,1 ),
                    x     =c(5,4,3,1,1,1,2,2,2,4,6,1,2 ))

# When creating the censoring time distribution remember that
#  censors happen after deaths, so the distribution does not drop until
#  time 3+, 4+, 5+, 6+ and 9+
csurv <- list(time=c(0, 3, 4, 5, 6, 9),
              p = cumprod(c(1, 10/11, 8/9, 5/7, 4/5, 1/2)))
#
# For estimation of event type 1, the first subject of event type
#  2 will have weights of curve$p over (0,3], (3,4], (4,5], (5,6], (6,9] 
#  and 9+.  All that really matters is the weight at times 1, 4, 5,
#  and 10, however, which are the points at which events of type 1 happen.
#  
# The next subjects of event type 2 occur at time 7 and 8, and will have a
#  weight of S(10)/S(6) = 1/2 at the time 10 event.
# For events of type 2 is the subjects of type 1 who get extended.

tid <- c(1, 2,2,2,2, 3:9, 10, 10, 11, 11, 12, 13)
temp1 <- data.frame(id=tid, start=0, fdata[tid,],
                    wt= c(1, csurv$p[c(1,2,3,6)], 1,1,1,1,1,1,1, 1, 
                         .5, 1, .5, 1,1),
                    endpoint=1)
temp1$time[c(2,3,4,5, 13,14, 15, 16)] <- c(3,4,5,10, 5,10, 5,10)
temp1$start[c(2,3,4,5, 13,14, 15,16)] <- c(0,3,4,5, 0,5, 0,5)
row.names(temp1) <- NULL
tfit <- survfit(Surv(start, time, status==1) ~1, temp1, weight=wt)

stat2 <- factor(fdata$status, 0:2, c("cen", "type1", "type2"))
temp2 <- finegray(Surv(time, stat2) ~x, fdata)


ctest <- crprep("time", "status", data=fdata, trans=1:2, cens=0, 
                id=1:13, keep='x')
xfit <- survfit(Surv(Tstart, Tstop, status==1) ~1, ctest, weight=weight.cens,
                subset=(failcode==1))

