# Tests of pseudovalues, by calculating directly from survfit and residuals
#  this assumes that residuals.survfit is correct
library(survival)
aeq <- function(x, y, ...) all.equal(as.vector(x), as.vector(y), ...)

mdata <- mgus2
temp <- ifelse(mdata$pstat==1, 1, 2*mdata$death)
mdata$event <- factor(temp, 0:2, c("censor", "pcm", "death"))
mdata$etime <- ifelse(mdata$pstat==1, mdata$ptime, mdata$futime)
mdata <- subset(mdata, etime > 12)  # remove first year

# Single endpoint, one curve
fit1 <- survfit(Surv(ptime, pstat) ~1, mdata)
# a time point before first event, after last event, at an event time,
#  and between event times
tvec <- c(10, 100, 250, 450)
rr1 <- resid(fit1, tvec)
sv1 <- summary(fit1, time=tvec, extend=TRUE)$surv

# one time point  
ps1a <- pseudo(fit1, time=100)
aeq(ps1a, sv1[2] + fit1$n*rr1[,2])
# multiple
ps1b <- pseudo(fit1,  time=c(10, 100, 250, 450))
aeq(ps1b,  sv1[col(rr1)] + fit1$n * rr1)

# Single endpoint, multiple curves
fit2 <- survfit(Surv(futime, death) ~ sex, mdata)
rr2 <- resid(fit2, time=tvec)
sv2 <- summary(fit2, time=tvec, extend=TRUE)$surv
sv2 <- t(matrix(sv2, ncol=2))   # row 1= female, row2 = male

# residuals are the same as for separate models
fit2a <- survfit(Surv(futime, death) ~1, mdata, subset=( sex=='F'))
fit2b <- survfit(Surv(futime, death) ~1, mdata, subset= (sex=='M'))
fem <- (mdata$sex=='F')
rr2a <- resid(fit2a, times=tvec)
rr2b <- resid(fit2b, times=tvec)
all.equal(rr2a, rr2[fem,])
all.equal(rr2b, rr2[!fem,])

# one time point
ps2a <- pseudo(fit2a, time=100)
aeq(ps2a, sv2[1,2] + fit2a$n[1]* rr2a[,2])
ps2b <- pseudo(fit2b, time=100)
aeq(ps2b, sv2[2,2] + fit2b$n[1]* rr2b[,2])

# overall psuedo are the same as for separate models
#  (each row of mdata belongs to a single curve)
ps2c <- pseudo(fit2, time=100)
aeq(ps2c[ fem], ps2a)
aeq(ps2c[!fem], ps2b)

# multiple time points
ps2d <- pseudo(fit2a, times=tvec)
aeq(ps2d, sv2[1, col(rr2a)] + fit2$n[1]* rr2a)
ps2e <- pseudo(fit2b, times=tvec)
aeq(ps2e, sv2[2, col(rr2b)] + fit2$n[2]* rr2b)

ps2f <- pseudo(fit2, times=tvec)
all.equal(ps2d, ps2f[ fem,])
all.equal(ps2e, ps2f[!fem,])

# Repeat the process for a multi-state model
fit3 <- survfit(Surv(etime, event) ~ sex, mdata)
fit3a <- survfit(Surv(etime, event) ~1, mdata, subset= (sex=='F'))
fit3b <- survfit(Surv(etime, event) ~1, mdata, subset= (sex=='M'))
rr3 <-  resid(fit3, times=tvec)
rr3a <- resid(fit3a, times=tvec)
rr3b <- resid(fit3b, times=tvec)
all.equal(rr3[fem,,], rr3a)
all.equal(rr3[!fem,,], rr3b)

ps3 <- pseudo(fit3, times=tvec)
ps3a <- pseudo(fit3a, times=tvec)
ps3b <- pseudo(fit3b, times=tvec)
all.equal(ps3[ fem,,], ps3a)
all.equal(ps3[!fem,,], ps3b)

sv3 <- summary(fit3, times=tvec, extend=TRUE)$pstate
sv3 <- array(sv3, dim=c(4,2,3))      #times, curve, state
# ps3a has dimensions (number obs in fit3a, 4 timepoints, 3 states)
#  to each of the 4x3 combinations we need to add the value of the
#  survival curve at that time.  A loop is easiest
temp1 <- array(0, dim= dim(rr3a))
temp2 <- array(0, dim= dim(rr3b))
for (i in 1:4) { # each of the 4 times
    for (j in 1:3) {  # each of the 3 endpoints
        temp1[, i,j] <- sv3[i,1,j] + fit3$n[1]*rr3a[,i,j]
        temp2[, i,j] <- sv3[i,2,j] + fit3$n[2]*rr3b[,i,j]
    }
}
aeq(temp1, ps3a)
aeq(temp2, ps3b)

###########################
# All again, just the same, for cumulative hazards
#  Though there are 2 of them, vs 3 states.
#
rr1 <- resid(fit1, tvec, type="cumhaz")
sv1 <- summary(fit1, time=tvec, extend=TRUE)$cumhaz

# one time point  
ps1a <- pseudo(fit1, time=100, type="cumhaz")
aeq(ps1a, sv1[2] + fit1$n*rr1[,2])
# multiple
ps1b <- pseudo(fit1,  time=c(10, 100, 250, 450), type="cumhaz")
aeq(ps1b,  sv1[col(rr1)] + fit1$n * rr1)

# Single endpoint, multiple curves
fit2 <- survfit(Surv(futime, death) ~ sex, mdata)
rr2 <- resid(fit2, time=tvec, type="cumhaz")
sv2 <- summary(fit2, time=tvec, extend=TRUE)$cumhaz
sv2 <- t(matrix(sv2, ncol=2))   # row 1= female, row2 = male

# residuals are the same as for separate models
rr2a <- resid(fit2a, times=tvec, type= "cumhaz")
rr2b <- resid(fit2b, times=tvec, type= "cumhaz")
all.equal(rr2a, rr2[fem,])
all.equal(rr2b, rr2[!fem,])

# one time point
ps2a <- pseudo(fit2a, time=100, type="cumhaz")
aeq(ps2a, sv2[1,2] + fit2a$n[1]* rr2a[,2])
ps2b <- pseudo(fit2b, time=100, type="cumhaz")
aeq(ps2b, sv2[2,2] + fit2b$n[1]* rr2b[,2])

# overall psuedo are the same as for separate models
#  (each row of mdata belongs to a single curve)
ps2c <- pseudo(fit2, time=100, type="cumhaz")
aeq(ps2c[ fem], ps2a)
aeq(ps2c[!fem], ps2b)

# multiple time points
ps2d <- pseudo(fit2a, times=tvec, type="cumhaz")
aeq(ps2d, sv2[1, col(rr2a)] + fit2$n[1]* rr2a)
ps2e <- pseudo(fit2b, times=tvec, type= "cumhaz")
aeq(ps2e, sv2[2, col(rr2b)] + fit2$n[2]* rr2b)

ps2f <- pseudo(fit2, times=tvec, type="cumhaz")
all.equal(ps2d, ps2f[ fem,])
all.equal(ps2e, ps2f[!fem,])

# Repeat the process for a multi-state model
rr3 <-  resid(fit3, times=tvec, type="cumhaz")
rr3a <- resid(fit3a, times=tvec, type="cumhaz")
rr3b <- resid(fit3b, times=tvec, type="cumhaz")
all.equal(rr3[fem,,], rr3a)
all.equal(rr3[!fem,,], rr3b)

ps3 <- pseudo(fit3, times=tvec, type="cumhaz")
ps3a <- pseudo(fit3a, times=tvec, type="cumhaz")
ps3b <- pseudo(fit3b, times=tvec, type="cumhaz")
all.equal(ps3[ fem,,], ps3a)
all.equal(ps3[!fem,,], ps3b)

sv3 <- summary(fit3, times=tvec, extend=TRUE)$cumhaz
sv3 <- array(sv3, dim=c(4,2,2))      #times, curve, state
# ps3a has dimensions (number obs in fit3a, 4 timepoints, 3 states)
#  to each of the 4x3 combinations we need to add the value of the
#  survival curve at that time.  A loop is easiest
temp1 <- array(0, dim= dim(rr3a))
temp2 <- array(0, dim= dim(rr3b))
for (i in 1:4) { # each of the 4 times
    for (j in 1:2) {  # each of the 3 endpoints
        temp1[, i,j] <- sv3[i,1,j] + fit3$n[1]*rr3a[,i,j]
        temp2[, i,j] <- sv3[i,2,j] + fit3$n[2]*rr3b[,i,j]
    }
}
aeq(temp1, ps3a)
aeq(temp2, ps3b)

#################################################
# Last, one more time with AUC
#  A bit more bother, since summary.survfit only returns AUC for one time
#   value at a time.  It also does not like times before the first event
#
tvec <- tvec[2:4]

rr1 <- resid(fit1, tvec, type="auc")
afun <- function(fit, times) {
    ntime <- length(times)
    if (length(fit$strata)) xfun <- function(x) x$table[, "*rmean"]
        else xfun <- function(x) x$table["*rmean"]

    temp <- xfun(summary(fit, rmean=times[1]))
    if (ntime==1) return(temp)
    
    result <- matrix(0, ntime, length(temp))
    result[1,] <- temp
    for (i in 2:ntime) 
        result[i,] <- xfun(summary(fit, rmean=times[i]))
    drop(result)
}
    
sv1 <- afun(fit1, tvec)

# one time point  
ps1a <- pseudo(fit1, time=100, type="auc")
aeq(ps1a, sv1[1] + fit1$n*rr1[,1])
# multiple
ps1b <- pseudo(fit1,  time=tvec, type="auc")
aeq(ps1b,  sv1[col(rr1)] + fit1$n * rr1)

# Single endpoint, multiple curves
rr2 <- resid(fit2, time=tvec, type="auc")
sv2 <- t(afun(fit2, tvec))

# residuals are the same as for separate models
rr2a <- resid(fit2a, times=tvec, type= "auc")
rr2b <- resid(fit2b, times=tvec, type= "auc")
all.equal(rr2a, rr2[fem,])
all.equal(rr2b, rr2[!fem,])

# one time point
ps2a <- pseudo(fit2a, time=100, type="auc")
aeq(ps2a, sv2[1,1] + fit2a$n[1]* rr2a[,1])
ps2b <- pseudo(fit2b, time=100, type="auc")
aeq(ps2b, sv2[2,1] + fit2b$n[1]* rr2b[,1])

# overall psuedo are the same as for separate models
#  (each row of mdata belongs to a single curve)
ps2c <- pseudo(fit2, time=100, type="auc")
aeq(ps2c[ fem], ps2a)
aeq(ps2c[!fem], ps2b)

# multiple time points
ps2d <- pseudo(fit2a, times=tvec, type="auc")
aeq(ps2d, sv2[1, col(rr2a)] + fit2$n[1]* rr2a)
ps2e <- pseudo(fit2b, times=tvec, type= "auc")
aeq(ps2e, sv2[2, col(rr2b)] + fit2$n[2]* rr2b)

ps2f <- pseudo(fit2, times=tvec, type="auc")
all.equal(ps2d, ps2f[ fem,])
all.equal(ps2e, ps2f[!fem,])

# Repeat the process for a multi-state model
rr3 <-  resid(fit3, times=tvec, type="auc")
rr3a <- resid(fit3a, times=tvec, type="auc")
rr3b <- resid(fit3b, times=tvec, type="auc")
all.equal(rr3[fem,,], rr3a)
all.equal(rr3[!fem,,], rr3b)

ps3 <- pseudo(fit3, times=tvec, type="auc")
ps3a <- pseudo(fit3a, times=tvec, type="auc")
ps3b <- pseudo(fit3b, times=tvec, type="auc")
all.equal(ps3[ fem,,], ps3a)
all.equal(ps3[!fem,,], ps3b)

sv3 <- rbind(summary(fit3, rmean=tvec[1])$table[,3],
             summary(fit3, rmean=tvec[2])$table[,3],
             summary(fit3, rmean=tvec[3])$table[,3])
sv3 <- array(sv3, dim=c(3,2,3))      #times, curve, state
# ps3a has dimensions (number obs in fit3a, 4 timepoints, 3 states)
#  to each of the 4x3 combinations we need to add the value of the
#  survival curve at that time.  A loop is easiest
temp1 <- array(0, dim= dim(rr3a))
temp2 <- array(0, dim= dim(rr3b))
for (i in 1:3) { # each of the 3 times
    for (j in 1:3) {  # each of the 3 states
        temp1[, i,j] <- sv3[i,1,j] + fit3$n[1]*rr3a[,i,j]
        temp2[, i,j] <- sv3[i,2,j] + fit3$n[2]*rr3b[,i,j]
    }
}
aeq(temp1, ps3a)
aeq(temp2, ps3b)
