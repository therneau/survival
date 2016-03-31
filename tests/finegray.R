library(survival)

# Test data set for Fine-Gray regression
fdata <- data.frame(time  =c(1,2,3,4,4,4,5,5,6,8,8, 9,10,12),
                    status=factor(c(1,2,0,1,0,0,2,1,0,0,2, 0,1 ,0), 0:2,
                             c("cen", "type1", "type2")),      
                    x     =c(5,4,3,1,2,1,1,2,2,4,6,1,2, 0),
                    id = 1:14)
test1 <- finegray(Surv(time, status) ~., fdata, count="fgcount")
test2 <- finegray(Surv(time, status) ~x, fdata, endpoint="type2")

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

all.equal(test1$id, c(1, 2,2,2,2, 3:6, 7, 7, 8:11, 11, 12:14))
twt <- c(1, csurv$p[c(1,2,3,6)], 1,1,1, 1, 1, 5/12, 1,1,1,  
                         1, 1/2, 1,1,1)
all.equal(test1$fgwt, twt)
#extra obs will end at times found in csurv$time, or max(time)=12
all.equal(test1$fgt2[test1$fgcount>0], c(4,6,12, 12,12))

#
# Verify the data reproduces a multi-state curve
#  censoring times may be different in the two setups so only 
#  compare at the event times
sfit <- survfit(Surv(time, status) ~1, fdata)
sfit1<- survfit(Surv(fgtime1, fgtime2, fgstatus) ~1, test1, weight=fgwt)
i1 <- sfit$n.event[,1] > 0
i2 <- sfit1$n.event > 0
all.equal(sfit$pstate[i1, 1], 1- sfit1$surv[i2])

sfit2 <- survfit(Surv(fgtime1, fgtime2, fgstatus) ~1, test2, weight=fgwt)
i1 <- sfit$n.event[,2] > 0
i2 <- sfit2$n.event > 0
all.equal(sfit$pstate[i1, 2], 1- sfit2$surv[i2])


