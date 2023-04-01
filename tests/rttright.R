library(survival)

# start with the example used in chapter 2 of the book

bdata <- data.frame(time =   c(1, 2, 2, 3, 4, 4, 5, 5, 8, 8, 
                               9, 10,11, 12,14, 15, 16, 16, 18, 20),
                    status = c(1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1,
                               0, 0, 1, 0, 0, 1, 0, 1, 0))

# First check: verify that the the RTTR reproduces the KM
kfit <- survfit(Surv(time, status) ~1, bdata)
bwt  <- rttright(Surv(time, status) ~1, bdata, renorm= FALSE)

cdf <- cumsum(bwt)/nrow(bdata)  # weighted CDF
cdf <- cdf[!duplicated(bdata$time, fromLast=TRUE)]  # remove duplicates
all.equal(kfit$surv, 1-cdf)


# A covariate divides both survfit and rttr into disjoint groups, so repeat
#  the above check on subsets of the aml data
afit <- survfit(Surv(time, status) ~x, aml)
awt <-  rttright(Surv(time, status) ~x, aml, renorm=TRUE)

igroup <- as.numeric(aml$x)
for (i in 1:2) {
    atemp <- awt[igroup ==i]   # subset for this curve
    index <- order(aml$time[igroup ==i])
    acdf <- cumsum(atemp[index])
    acdf <- acdf[!duplicated(aml$time[igroup ==i], fromLast=TRUE)]
    print(all.equal(afit[i]$surv, 1-acdf))
}


###########
# Alternate computation using inverse prob of censoring weights.
# First shift the censorings to avoid ties: if there is a death and a censor
#   at time 10, say, the death was not at risk of censoring. Censoring weights
#   happen "later".  This also results in a left-continuous curve.
delta <- min(diff(sort(unique(bdata$time)))) /3
offset <- ifelse(bdata$status==1, 0, delta)
cfit <- survfit(Surv(time+ offset, 1-status) ~ 1, bdata)

# interpolate
indx <- findInterval(bdata$time, cfit$time)
cwt <- ifelse(bdata$status==0, 0, 1/cfit$surv[indx])
all.equal(bwt, cwt) 

# Multiple time points, this example is used in the vignette
tdata <- data.frame(time=   c(1,2,2,3,4,4,5,5,8,9),
                    status= c(1,1,0,1,0,0,1,0,1,1))
fit1 <- rttright(Surv(time, status) ~ 1, tdata, times=2:6, renorm=FALSE)
fit2 <- rttright(Surv(time, status) ~ 1, tdata, times=2:6, renorm=TRUE)
all.equal(fit1, 10*fit2)
all.equal(fit1, cbind(7, c(7,7,0,8,8,8,8,8,8,8), 
                         c(7,7,0,8,8,8,8,8,8,8), 
                         c(7,7,0,8,0,0,12,12,12,12),
                         c(7,7,0,8,0,0,12, 0, 18,18))/7, check.attributes=FALSE)

# Now test with (start, stop] data, should get the same results
b2 <- survSplit(Surv(time, status) ~ 1, bdata, cut= c(3,5, 7, 14),
                id = "subject")
indx <- c(seq(1, 65, by=2), seq(64, 2, by= -2))
b2 <- b2[indx,]    # not in time within subject order (stronger test)

b2wt <- rttright(Surv(tstart, time, status) ~1, b2, id=subject)
indx2 <- order(b2$time)
cdf2 <- cumsum(b2wt[indx2])
cdf2 <- cdf2[!duplicated(b2$time[indx2], fromLast=TRUE)] # remove duplicates
utime2 <- sort(unique(b2$time))   # will have an extra time 7
utime1 <- sort(unique(bdata$time))
all.equal(cdf2[match(utime1, utime2)], cdf)


# Competing risks
mdata <- mgus2
mdata$etime <- with(mgus2, ifelse(pstat==1, ptime, futime))
mdata$estat <- with(mgus2, ifelse(pstat==1, 1, 2*death))
mdata$estat <- factor(mdata$estat, 0:2, c('censor', 'pcm', 'death'))
mfit <- survfit(Surv(etime, estat) ~1, mdata, id=id)
mwt1 <- rttright(Surv(etime, estat) ~1, mdata, id=id)

morder <- order(mdata$etime)
mdata2 <- mdata[morder,]
mwt2   <- rttright(Surv(etime,estat) ~1, mdata2, id=id)
all.equal(mwt1[morder], mwt2)

keep <- !duplicated(mdata2$etime, fromLast=TRUE)
csum1 <- cumsum(ifelse(mdata2$estat=="pcm", mwt2, 0))
csum2 <- cumsum(ifelse(mdata2$estat=="death", mwt2, 0))

all.equal(mfit$pstate[,2], csum1[keep])
all.equal(mfit$pstate[,3], csum2[keep])

# Case weights, at multiple times
bwt <- rep(1:2, length=nrow(bdata))
tm <- c(2, 6, 10, 15, 18)
fit1 <- rttright(Surv(time, status) ~1, bdata, weights=bwt, times= tm)
casefit <- survfit(Surv(time, status) ~ 1, bdata, weights= bwt)
csum1 <- summary(casefit, censor=FALSE, times= tm)
for (i in 1:length(tm)) {
    c1 <- sum(fit1[bdata$status==1 & bdata$time <= tm[i], i])
    print(all.equal(c1, 1-csum1$surv[i]))
}


###
# Delayed entry, tiny data set
#  Still to be done
delay <- data.frame(t0=c(0,0,0,0,3,0),
                    t1=1:6,
                    status=c(1,0,1,0,0,1),
                    id=1:6)
# dwt <- rttright(Surv(t0, t1, status) ~ 1, delay, id=id, times=0:5 + .9)


