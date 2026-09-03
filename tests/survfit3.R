# Check that three version: Surv2, Surv, and status as a factor all give the
#  same answer.  These are different paths through the code, and the last
#  uses survfitaj.c rather than survfitkm.c
library(survival)
aeq <- function(x, y, ...) all.equal(as.vector(x), as.vector(y), ...)

tdata <- data.frame(id = c(1, 1,  1,  1, 2, 2, 3, 3, 4, 4,  4, 5, 5, 5, 5),
                    tt = c(0, 5, 10, 15, 4, 6, 1, 4, 5, 9, 11, 0, 3, 4, 10),
                    ss=  c(0, 0,  0,  1, 0, 0, 0, 0, 0, 0, 1,  0, 0, 0, 1),
                    x  = c(1, 1,  1,  1, 2, 2, 1, 2, 1, 1, 2,  1, 2, 2, 20))
cdata <- fromtimeline(Surv2(tt,ss) ~ ., tdata, id=id)

s1 <- survfit(Surv2(tt, ss) ~ x, tdata, id=id)
s2 <- survfit(Surv(tt1, tt2, ss) ~x, cdata, id=id)
s3 <- survfit(Surv(tt1, tt2, factor(ss==1)) ~ x, cdata, id=id)
aeq(s1$surv, s2$surv)
aeq(s1$surv, s3$pstate[,1])
aeq(s1$time, c(3,4,11,15, 6, 10))

# The timelines are (0, 5, 10, 15), (1,4), (5,9,11) and (0,3) for x=1,
#    and (4,6) (3,4,10) for x=2, all simple survival.
# Subject 5 has a foot in each camp (x=0 and x=1), 3 and 4 appear to but
#   actually don't as there is no further follow-up after x changes.
#
# The in-between times are not censoring points, so the time points for s1,
#  s2, s3 should be (3, 4, 11, and 15) for x=1, and (6, 10) for x=2. This
#  test arose from a case where the survfitKM path (s1 and s2) was including
#  some of them but the survfitAJ path was not; found in a much larger data set.
#
# With entry=TRUE include the entry points of 0, 1, 5 for x=1 and 4,6 for x=2
s4 <- survfit(Surv2(tt, ss) ~ x, tdata, id=id, entry=TRUE)
s5 <- survfit(Surv(tt1, tt2,ss) ~ x, cdata, id=id, entry=TRUE)
s6 <- survfit(Surv(tt1, tt2, factor(ss==1)) ~ x, cdata, id=id, entry=TRUE)
aeq(s4$time, c(0,1,3,4,5,11,15, 3,4,6,10))
aeq(s4$time, s5$time)
aeq(s6$time, s5$time)
aeq(s4$surv[s4$n.enter ==0], s1$surv)
aeq(s4$surv, s5$surv)
aeq(s4$surv, s6$pstate[,1])

