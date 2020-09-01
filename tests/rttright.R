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
    cdf <- cumsum(atemp[index])/length(atemp)
    cdf <- cdf[!duplicated(aml$time[igroup ==i], fromLast=TRUE)]
    print(all.equal(afit[i]$surv, 1-cdf))
}

