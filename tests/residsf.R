library(survival)
aeq <- function(x, y, ...) all.equal(as.vector(x), as.vector(y), ...)
#
# Tests of the residuals.survfit function
#
# The influence argument of survfit returns all the residuals at every time
#  point, but for large data sets the result will be huge.  This function uses
#  a different algorithm which will be faster when the number of time
#  points being reported out is small. 

# Start with small data sets and work up.  First simple survival.
test1 <- data.frame(time=  c(9, 3,1,1,6,6,8),
                    status=c(1,NA,1,0,1,1,0),
                    x=     c(0, 2,1,1,1,0,0))
indx <- order(test1$time[!is.na(test1$status)])

s1 <- survfit(Surv(time, status) ~1, test1, influence=3)
# true influence for survival and hazard, in time order
inf1 <- matrix(c(-20, rep(4,5), -10, 2, -13, -13, 17, 17,
                 rep(0,6))/144, ncol=3,
               dimnames=list(1:6, c(1,6,9)))
inf2 <- matrix(c(10, rep(-2,5), 10, -2, 7,7, -11, -11)/72,
               ncol=2)
aeq(s1$influence.surv[indx,], inf1[, c(1,2,2,3)])
aeq(s1$influence.chaz[indx,], inf2[,c(1,2,2,2)])

r1 <- resid(s1, times=c(0, 3, 5, 8, 10))
all(r1[,1] ==0)
aeq(r1[indx,2:5], inf1[,c(1,1,2,3)])

r2 <- resid(s1, times=c(0, 3, 5, 8, 10), type="cumhaz")
all(r2[,1] ==0)
aeq(r2[indx,2:5], inf2[,c(1,1,2,2)])

# AUC is a sum of rectangles, height= S, width based on time points,
#  so the leverage is a weighted sum of dfbeta values for S
r3 <- resid(s1, times=c(1,4, 8, 10), type="sojourn")
inf3 <- inf1 %*% cbind(c(0,0,0), c(3,0,0), c(5,2,0), c(5,3,1))
aeq(r3[indx,], inf3)

# exp(Nelson-Aalen) 
s2 <- survfit(Surv(time, status) ~1, test1, stype=2, influence=3)
r4 <- resid(s2, times=c(0, 3, 5, 8, 10), type="pstate")
inf4 <- -inf2[, c(1,2,2)] %*% diag(s2$surv[c(1,2,4)])
aeq(r4[indx,2:5], inf4[,c(1,1,2,3)])
aeq(s2$influence.surv[indx,], inf4[,c(1,2,2,3)])

r5 <- resid(s2, times=c(1,4, 8, 10), type="sojourn") 
inf5 <- inf4 %*% cbind(c(0,0,0), c(3,0,0), c(5,2,0), c(5,3,1))
aeq(r5[indx,], inf5)

# Fleming-Harrington 
# This one is hard, the code still fails
s3 <- survfit(Surv(time, status) ~1, test1, ctype=2, influence=2)
inf6 <-  matrix(c( rep(c(5, -1), c(1, 5))/36, c(5,-1)/36, 
                  c(21,21,-29, -29)/144), ncol=2)
# r6 <- resid(s3, times =c(0, 3, 5, 8, 10), type="cumhaz")

# Part 2: single state, with start/stop data, multiple curves,
#   second curve is identical to test1
#   Then put it out of order.  

test2 <- data.frame(t1  =c(1, 2, 5, 2, 1, 7, 3, 4, 8, 8,
                           0,0,0,0,0,0),
                    t2  =c(2, 3, 6, 7, 8, 9, 9, 9,14, 17,
                           9, 1, 1, 6, 6, 8),
                    event=c(1, 1, 1, 1, 1, 1, 1, 0, 0, 0,
                            1, 1, 0, 1, 1, 0),
                    x    = rep(1:2, c(10, 6)),
                    id   = 1:16)

s4 <- survfit(Surv(t1, t2, event) ~ x, test2, influence=TRUE)
r6 <- resid(s4, time=c(4, 8, 10), type="surv")
aeq(r6[1:10,], s4$influence.surv[[1]][,c(2, 5, 6)])
aeq(r6[11:16,],s4$influence.surv[[2]][,c(1,3, 4)])
aeq(r6[11:16,2:3], r1[,4:5])

r7 <- resid(s4, time=c(4, 8, 10), type="cumhaz")
aeq(r7[1:10,], s4$influence.chaz[[1]][,c(2, 5, 6)])
aeq(r7[11:16,],s4$influence.chaz[[2]][,c(1,3, 4)])
aeq(r7[11:16, 2:3], r2[,4:5])

# Compute the AUC at times 8 and 10, the first is a reporting time, the
#  second is in between
r8 <- resid(s4, time= c(8, 10), type="auc")
aeq(r8[11:16,], r3[,3:4])

# curve1:
inf1 <- s4$influence.surv[[1]]
d1 <- inf1[,1:4] %*% diff(s4$time[1:5])
d2 <- inf1[,1:6] %*% diff(c(s4$time[1:6], 10))
aeq(cbind(d1, d2), r8[1:10,])

# curve2:
inf2 <- s4$influence.surv[[2]]
d3 <- inf2[,1:2] %*% diff(s4$time[9:11])
d4 <- inf2[,1:4] %*% diff(c(s4$time[9:12], 10))
aeq(cbind(d3, d4), r8[11:16,])

# scramble the data
reord <- c(1,3,5,7,9,11,13, 15,2,4,6,8,10,12,14,16)
test2b <-test2[reord,]
s5 <- survfit(Surv(t1, t2, event) ~x, test2b, influence=TRUE)
r9 <- resid(s5, time=c(4, 8, 10), type="surv")
aeq(r6[reord,], r9)
 
# 
# For multistate use the same data set as mstate.R, where results have been
#  worked out by hand.
#
tdata <- data.frame(id= LETTERS[3*c(1, 1, 1,  2,  3,  4, 4, 4,  5,  5)],
                    t1= c(0, 4, 9,  1,  2,  0, 2, 8,  1,  3),
                    t2= c(4, 9, 10, 5,  9,  2, 8, 9,  3, 11),
                    st= c(1, 2,  1, 2,  3,  1, 3, 0,  3,  0),
                    i0= c(1, 2,  3, 2,  1,  1, 2, 4,  3,  4),
                    wt= 1:10)

tdata$st <- factor(tdata$st, c(0:3),
                    labels=c("censor", "a", "b", "c"))
tdata$i0 <- factor(tdata$i0, 1:4,
                    labels=c("entry", "a", "b", "c"))  
 
tfun <- function(data=tdata) {
    reorder <- c(10, 9, 1, 2, 5, 4, 3, 7, 8, 6)
    new <- data[reorder,]
    new
}
mtest2 <- tfun(tdata)  # scrambled version

mfit1 <- survfit(Surv(t1, t2, st) ~ 1, tdata, id=id, istate=i0,
                 influence=1)

test1a <- resid(mfit1, time=c(3, 7, 9), method=1)
#test1b <- resid(mfit1, time=c(3, 7, 9), method=2)
#aeq(test1a, test1b)
aeq(test1a, mfit1$influence.pstate[,c(3,5,7),])

# AUC, start simple - auc at final time
test3 <- resid(mfit1, time=11, type='RMST')
delta <- diff(c(0, mfit1$time))
s1 <- apply(mfit1$influence[, 1:8, ], c(1,3), function(x) sum(delta*x))
aeq(test3, s1)

# extend to an earlier and later time
test3b <- resid(mfit1, time=c(-1,11,15), type='RMST')
all(test3b[,1,] ==0)
aeq(test3b[,2,], s1)
aeq(test3b[,3,], s1 + mfit1$influence[,9,]*4)

#
# competing risks
#
mdata <- mgus2
mdata$etime <- with(mdata, ifelse(pstat==1, ptime, futime))
temp <- with(mdata, ifelse(pstat==1, 1, 2*death))
mdata$event <- factor(temp, 0:2, c("censor", "PCM", "Death"))
mfit <- survfit(Surv(etime, event) ~1, mdata, influence=1)
rr <- resid(mfit, time=360)
index <- sum(mfit$time <= 360) +1  # influence has time 0 too
aeq(mfit$influence.pstate[,index,], rr)

