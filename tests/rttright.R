library(survival)

# start with the example used in chapter 2 of the book

bdata <- data.frame(time =   c(1, 2, 2, 3, 4, 4, 5, 5, 8, 8, 
                               9, 10,11, 12,14, 15, 16, 16, 18, 20),
                    status = c(1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1,
                               0, 0, 1, 0, 0, 1, 0, 1, 0))

kfit <- survfit(Surv(time, status) ~1, bdata)
bwt  <- rttright(Surv(time, status) ~1, bdata)

cdf <- cumsum(bwt)/nrow(bdata)  # weighted CDF
cdf <- cdf[!duplicated(bdata$time, fromLast=TRUE)]  # remove duplicates
all.equal(kfit$surv, 1-cdf)


afit <- survfit(Surv(time, status) ~x, aml)
awt <-  rttright(Surv(time, status) ~x, aml)

igroup <- as.numeric(aml$x)
for (i in 1:2) {
    atemp <- awt[igroup ==i]   # subset for this curve
    index <- order(aml$time[igroup ==i])
    acdf <- cumsum(atemp[index])/length(atemp)
    acdf <- acdf[!duplicated(aml$time[igroup ==i], fromLast=TRUE)]
    print(all.equal(afit[i]$surv, 1-acdf))
}

# Now test with (start, stop] data
b2 <- survSplit(Surv(time, status) ~ 1, bdata, cut= c(3,5, 7, 14),
                id = "subject")
indx <- c(seq(1, 65, by=2), seq(64, 2, by= -2))
b2 <- b2[indx,]    # not in time within subject order

b2wt <- rttright(Surv(tstart, time, status) ~1, b2, id=subject)
indx2 <- order(b2$time)
cdf2 <- cumsum(b2wt[indx2])/length(unique(b2$subject))
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
csum1 <- cumsum(ifelse(mdata2$estat=="pcm", mwt2, 0))/nrow(mdata2)
csum2 <- cumsum(ifelse(mdata2$estat=="death", mwt2, 0))/nrow(mdata2)

all.equal(mfit$pstate[,2], csum1[keep])
all.equal(mfit$pstate[,3], csum2[keep])

###
# Delayed entry, tiny data set
#  This is the data that showed me that the RTTR idea does not extend
#  to delayed entry, compute it how you will.
delay <- data.frame(t0=c(0,0,0,0,3,0),
                    t1=1:6,
                    status=c(1,0,1,0,0,1),
                    id=1:6)
# dwt <- rttright(Surv(t0, t1, status) ~ 1, delay, id=id,
#                 times=0:5 + .9)

