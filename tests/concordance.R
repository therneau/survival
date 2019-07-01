library(survival)
options(na.action=na.exclude) # preserve missings
options(contrasts=c('contr.treatment', 'contr.poly')) #ensure constrast type

#
# Simple tests of survConcordance, which has been replaced by the more general
#  concordance function.
#
aeq <- function(x,y, ...) all.equal(as.vector(x), as.vector(y), ...)

grank <- function(x, time, grp, wt) 
    unlist(tapply(x, grp, rank))
grank2 <- function(x, time, grp, wt) {  #for case weights
    if (length(wt)==0) wt <- rep(1, length(x))
    z <- double(length(x))
    for (i in unique(grp)) {
        indx <- which(grp==i)
        temp <- tapply(wt[indx], x[indx], sum)
        temp <- temp/2  + c(0, cumsum(temp)[-length(temp)])
        z[indx] <- temp[match(x[indx], names(temp))]
    }
    z
}

# Concordance by brute force.  O(n^2) algorithm, but ok for n<500 or so
allpair <- function(x, time, status, wt, all=FALSE) {
    if (missing(wt)) wt <- rep(1, length(x))
    count <- sapply(which(status==1), function(i) {
        atrisk <- (time > time[i]) | (time==time[i] & status==0)
        temp <- tapply(wt[atrisk], factor(sign(x[i] -x[atrisk]), c(1, -1, 0)),
                       sum)
        wt[i]* c(ifelse(is.na(temp), 0, temp),
                 (sum(wt[time==time[i] & status==1]) - wt[i])/2) 
    })
    rownames(count) <- c("concordant", "discordant", "tied.x", "tied.y")
    if (all) {
        colnames(count) <- time[status==1]
        t(count)
    }
    else rowSums(count)
}


# The std of C = std(numerator)/(number of comparable pairs)
# The information matrix of a Cox model is = to the var(C-D)
cfun <- function(fit) fit$cvar * sum(fit$count[1:3])^2
                
tdata <- aml[aml$x=='Maintained', c("time", "status")]
tdata$x <- c(1,6,2,7,3,7,3,8,4,4,5)
tdata$wt <- c(1,2,3,2,1,2,3,4,3,2,1)
fit <- concordance(Surv(time, status) ~x, tdata)

aeq(fit$count[1:4], c(24,14,2,0))
cfit <- coxph(Surv(time, status) ~ tt(x), tdata, tt=grank, method='breslow',
              iter=0, x=T)
cdt <- coxph.detail(cfit)
aeq(sum(cdt$imat), cfun(fit))
aeq(sum(2*cdt$score), diff(fit$count[1:2]))
aeq(with(tdata, allpair(x, time, status)), c(14,24,2,0))

# Lots of ties
tempy <- Surv(c(1,2,2,2,3,4,4,4,5,2), c(1,0,1,0,1,0,1,1,0,1))
tempx <- c(5,5,4,4,3,3,7,6,5,4)
fit2 <- concordance(tempy ~ tempx)
addxy <- function(x) c(x[1:3], sum(x[4:5]))
aeq(addxy(fit2$count), allpair(tempx, tempy[,1], tempy[,2]))
cfit2 <-  coxph(tempy ~ tt(tempx), tt=grank, method='breslow', iter=0)
aeq(cfit2$var, 1/cfun(fit2))

# Bigger data
fit3 <- concordance(Surv(time, status) ~ age, lung, reverse=TRUE)
aeq(addxy(fit3$count), allpair(lung$age, lung$time, lung$status-1))
cfit3 <- coxph(Surv(time, status) ~ tt(age), lung, 
               iter=0, method='breslow', tt=grank, x=T)
cdt <- coxph.detail(cfit3)
aeq(sum(cdt$imat), cfun(fit3)) 
aeq(2*sum(cdt$score), diff(fit3$count[2:1]))


# More ties
fit4 <- concordance(Surv(time, status) ~ ph.ecog, lung, reverse=TRUE)
aeq(addxy(fit4$count), allpair(lung$ph.ecog, lung$time, lung$status-1))
aeq(fit4$count[1:5], c(8392, 4258, 7137, 21, 7))
cfit4 <- coxph(Surv(time, status) ~ tt(ph.ecog), lung, 
               iter=0, method='breslow', tt=grank)
aeq(1/cfit4$var, cfun(fit4))

# Case weights
fit5 <- concordance(Surv(time, status) ~ x, tdata, weight=wt, reverse=TRUE)
fit6 <- concordance(Surv(time, status) ~x, tdata[rep(1:11,tdata$wt),])
aeq(addxy(fit5$count), with(tdata, allpair(x, time, status, wt)))
aeq(fit5$count[1:4], c(70, 91, 7, 0))  # checked by hand
aeq(fit5$count[1:3], fit6$count[c(2,1,3)])  #spurious "tied on time" values, ignore
aeq(fit5$std, fit6$std)
cfit5 <- coxph(Surv(time, status) ~ tt(x), tdata, weight=wt, 
               iter=0, method='breslow', tt=grank2)
cfit6 <- coxph(Surv(time, status) ~ tt(x), tdata[rep(1:11,tdata$wt),], 
               iter=0, method='breslow', tt=grank)
aeq(1/cfit6$var, cfun(fit6))
aeq(cfit5$var, cfit6$var)

# Start, stop simplest cases
fit7 <- concordance(Surv(rep(0,11), time, status) ~ x, tdata)
aeq(fit7$count, fit$count)
aeq(fit7$std.err, fit$std.err)
fit7 <- concordance(Surv(rep(0,11), time, status) ~ x, tdata, weight=wt)
aeq(fit5$count, fit7$count[c(2,1,3:5)])  #one reversed, one not

# Multiple intervals for some, but same risk sets as tdata
tdata2 <- data.frame(time1=c(0,3, 5,  6,7,   0,  4,17,  7,  0,16,  2,  0, 
                             0,9, 5),
                     time2=c(3,9, 13, 7,13, 18, 17,23, 28, 16,31, 34, 45, 
                             9,48, 60),
                     status=c(0,1, 1, 0,0,  1,  0,1, 0, 0,1, 1, 0, 0,1, 0),
                     x = c(1,1, 6, 2,2, 7, 3,3, 7, 3,3, 8, 4, 4,4, 5),
                     wt= c(1,1, 2, 3,3, 2, 1,1, 2, 3,3, 4, 3, 2,2, 1))
fit8 <- concordance(Surv(time1, time2, status) ~x, tdata2, weight=wt,
                    reverse=TRUE)
aeq(fit5$count, fit8$count)
aeq(fit5$std.err, fit8$std.err)
cfit8 <- coxph(Surv(time1, time2, status) ~ tt(x), tdata2, weight=wt, 
               iter=0, method='breslow', tt=grank2)
aeq(1/cfit8$var, cfun(fit8))

# Stratified
tdata3 <- data.frame(time1=c(tdata2$time1, rep(0, nrow(lung))),
                     time2=c(tdata2$time2, lung$time),
                     status = c(tdata2$status, lung$status -1),
                     x = c(tdata2$x, lung$ph.ecog),
                     wt= c(tdata2$wt, rep(1, nrow(lung))),
                     grp=rep(1:2, c(nrow(tdata2), nrow(lung))))
fit9 <- concordance(Surv(time1, time2, status) ~x + strata(grp),
                        data=tdata3, weight=wt, reverse=TRUE)
aeq(fit9$count[1,], fit5$count)
aeq(fit9$count[2,], fit4$count)
