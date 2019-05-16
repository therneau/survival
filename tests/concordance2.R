library(survival)
options(na.action=na.exclude) # preserve missings
options(contrasts=c('contr.treatment', 'contr.poly')) #ensure constrast type

#
# Tests for the condordance function. 
#
aeq <- function(x,y, ...) all.equal(as.vector(x), as.vector(y), ...)

grank <- function(x, time, grp, wt) {
    if (all(wt==1)) unlist(tapply(x, grp, rank))
    else unlist(tapply(1:length(x), grp, function(i) {
        xx <- x[i]  # x and wts for this subset of the data
        ww <- wt[i]
        temp <- outer(xx, xx, function(a, b) sign(b-a))
        colSums(ww*temp)/2
    }))
}

# a Cox model using iter=0, ties='breslow' and the above function has a score 
#  statistic which is U=(C-D)/2 and score test U^2/H, where H is the Cox model 
#  information matrix, with fit$var=1/H.  The concordance is U+ 1/2. 
# Pull out the Somers' d and its variance
phget <- function(fit) {
    c(d =  2*sqrt(fit$score/fit$var), v= 4/fit$var)
}
fscale <- function(fit) {
    if (is.matrix(fit$count)) temp <- colSums(fit$count) else temp <- fit$count
    npair <- sum(temp[1:3])
    c(d = abs(temp[1]-temp[2]), v=4*fit$cvar*npair^2)
}

# Concordance by brute force.  O(n^2) algorithm, but ok for n<500 or so
allpair <- function(time, status, x, wt, all=FALSE) {
    n <- length(time)
    if (missing(wt)) wt <- rep(1, length(x))
    count <- sapply(which(status==1), function(i) {
        atrisk <- (time > time[i]) | (time==time[i] & status==0)
        temp <- tapply(wt[atrisk], factor(sign(x[i] -x[atrisk]), c(-1, 1, 0)),
                       sum)
        tiedtime <- (time==time[i] & status ==1 & (1:n)>i)
        ties <- tapply(wt[tiedtime], factor(x[tiedtime]==x[i], 
                                            c(FALSE, TRUE)),sum)
        wt[i]* c(ifelse(is.na(temp), 0, temp), ifelse(is.na(ties), 0, ties))
    })
    rownames(count) <- c("concordant", "discordant", "tied.x", "tied.y",
                         "tied.xy")
    if (all) {
        colnames(count) <- time[status==1]
        t(count)
    }
    else rowSums(count)
}

# leverage by brute force
leverage <- function(time, status, x, wt, eps=1e-5) {
    if (missing(wt)) wt <- rep(1, length(x))
    toss <- is.na(time + status + x +wt)
    if (any(toss)) {
        time <- time[!toss]
        status <- status[!toss]
        x <- x[!toss]
        wt <- wt[!toss]
    }
    n <- length(time)
    influence <- matrix(0, n, 5)
    t2 <- time + eps*(status==0)
    for (i in 1:n) {
        if (status[i] ==0) comparable <- (time<=time[i] & status==1)
        else comparable <- ifelse(status==0, time >= time[i], time!=time[i])
        temp <- sign((x[i]-x[comparable])*(t2[i] - t2[comparable]))
        influence[i,1:3] <-tapply(wt[comparable],factor(temp, c(1,-1,0)), sum)
        if (status[i]==1) {
            tied <- (time==time[i] & status==1 & (1:n)!= i)
            if (any(tied)) {
                itemp<- tapply(wt[tied], factor(x[tied]==x[i], 
                                                c(FALSE, TRUE)), sum)
                influence[i,4:5] <- itemp
            }
        }
    }
    dimnames(influence) <- list(as.character(Surv(time, status)), 
                                c("concord", "discord", "tie.x", "tie.y",
                                  "tie.xy"))
    ifelse(is.na(influence), 0, influence)
}

tdata <- aml[aml$x=='Maintained', c("time", "status")]
tdata$x <- c(1,6,2,7,3,7,3,8,4,4,5)
tdata$wt <- c(1,2,3,2,1,2,3,4,3,2,1)

fit <- concordance(Surv(time, status) ~x, tdata, influence=2)
aeq(fit$count, with(tdata, allpair(time, status, x)))
aeq(fit$influence, with(tdata, leverage(time, status, x)))

cfit <- coxph(Surv(time, status) ~ tt(x), tdata, tt=grank, ties='breslow',
              iter=0, x=T)
aeq(phget(cfit), fscale(fit))  # agree with Cox model

# Test 2: Lots of ties
tempy <- Surv(c(1,2,2,2,3,4,4,4,5,2), c(1,0,1,0,1,0,1,1,0,1))
tempx <- c(5,5,4,4,3,3,7,6,5,4)
fit2 <- concordance(tempy ~ tempx, influence=2)
aeq(fit2$count, allpair(tempy[,1], tempy[,2], tempx))
aeq(fit2$influence, leverage(tempy[,1], tempy[,2], tempx))
cfit2 <- coxph(tempy ~ tt(tempx), tt=grank, ties="breslow", iter=0)
aeq(phget(cfit2), fscale(fit2))  # agree with Cox model

# Bigger data
cox3 <- coxph(Surv(time, status) ~ age + sex + ph.ecog, lung)
fit3 <- concordance(Surv(time, status) ~ predict(cox3), lung, influence=2)
aeq(fit3$count, allpair(lung$time, lung$status-1,predict(cox3)))
aeq(fit3$influence, leverage(lung$time, lung$status-1,predict(cox3)))
cfit3 <- coxph(Surv(time, status) ~ tt(predict(cox3)), tt=grank,
               ties="breslow", iter=0, data=lung)
aeq(phget(cfit3), fscale(fit3))  # agree with Cox model

# More ties
fit4  <- concordance(Surv(time, status) ~ ph.ecog, lung, influence=2)
fit4b <- concordance(Surv(time, status) ~ ph.ecog, lung, reverse=TRUE)
aeq(fit4$count, allpair(lung$time, lung$status-1, lung$ph.ecog))
aeq(fit4b$count, c(8392, 4258, 7137, 21, 7))
cfit4 <- coxph(Surv(time, status) ~ tt(ph.ecog), lung, 
               iter=0, method='breslow', tt=grank)
aeq(phget(cfit4), fscale(fit4))  # agree with Cox model

# Case weights
fit5 <- concordance(Surv(time, status) ~ x, tdata, weight=wt, influence=2)
fit6 <- concordance(Surv(time, status) ~x, tdata[rep(1:11,tdata$wt),])
aeq(fit5$count, with(tdata, allpair(time, status, x, wt)))
aeq(fit5$count, c(91, 70, 7, 0, 0))  # checked by hand
aeq(fit5$count[1:3], fit6$count[1:3])  #spurious "tied.xy" values, ignore
aeq(fit5$var[2], fit6$var[2])
aeq(fit5$influence, with(tdata, leverage(time, status, x, wt)))
cfit5 <- coxph(Surv(time, status) ~ tt(x), tdata, weight=wt, 
               iter=0, method='breslow', tt=grank)
aeq(phget(cfit5), fscale(fit5))  # agree with Cox model

# Start, stop simplest cases
fit6 <- concordance(Surv(rep(0,11), time, status) ~ x, tdata)
aeq(fit6$count, fit$count)
aeq(fit6$var, fit$var)
fit7 <- concordance(Surv(rep(0,11), time, status) ~ x, tdata, weight=wt)
aeq(fit7$count, fit5$count)
aeq(fit7$var, fit5$var)

# Multiple intervals for some, but same risk sets as tdata
tdata2 <- data.frame(time1=c(0,3, 5,  6,7,   0,  4,17,  7,  0,16,  2,  0, 
                             0,9, 5),
                     time2=c(3,9, 13, 7,13, 18, 17,23, 28, 16,31, 34, 45, 
                             9,48, 60),
                     status=c(0,1, 1, 0,0,  1,  0,1, 0, 0,1, 1, 0, 0,1, 0),
                     x = c(1,1, 6, 2,2, 7, 3,3, 7, 3,3, 8, 4, 4,4, 5),
                     wt= c(1,1, 2, 3,3, 2, 1,1, 2, 3,3, 4, 3, 2,2, 1),
                     id= c(1,1, 2, 3,3, 4, 5,5, 6, 7,7, 8, 9, 10,10, 11))
fit8 <- concordance(Surv(time1, time2, status) ~x, cluster=id, tdata2, 
                    weight=wt, influence=2)
aeq(fit5$count, fit8$count)
# influence has one row per obs, so the next line is false: mismatched lengths
# aeq(fit5$influence, fit8$influence) 
aeq(fit5$var, fit8$var)
cfit8 <- coxph(Surv(time1, time2, status) ~ tt(x), tdata2, weight=wt, 
               iter=0, method='breslow', tt=grank)
aeq(phget(cfit8), fscale(fit8))  # agree with Cox model

# Stratified
tdata3 <- data.frame(time1=c(tdata2$time1, rep(0, nrow(lung))),
                     time2=c(tdata2$time2, lung$time),
                     status = c(tdata2$status, lung$status -1),
                     x = c(tdata2$x, lung$ph.ecog),
                     wt= c(tdata2$wt, rep(1, nrow(lung))),
                     grp=rep(1:2, c(nrow(tdata2), nrow(lung))),
                     id = c(tdata2$id, 100+ 1:nrow(lung)))
fit9 <- concordance(Surv(time1, time2, status) ~x + strata(grp), cluster=id,
                        data=tdata3, weight=wt, influence=2)
aeq(fit9$count, rbind(fit8$count, fit4$count))

# check out case weights, strata, and grouped jackknife;
#   force several ties in x, y, and xy (with missing values too for good measure).
tdata <- subset(lung, select=-c(meal.cal, wt.loss, sex, age))
tdata$wt <- rep(1:25, length=nrow(tdata))/10
tdata$time <- ceiling(tdata$time/30)  # force ties in y
tfit <- coxph(Surv(time, status) ~ ph.ecog + pat.karno + strata(inst)
              + cluster(inst), tdata, weight=wt)
tdata$tpred <- predict(tfit)
cm4 <- concordance(tfit, influence=3, keepstrata=TRUE)
cm5 <- concordance(Surv(time, status) ~ tpred + strata(inst) + cluster(inst),
                   data=tdata, weight=wt, reverse=TRUE, influence=3,
                   keepstrata=TRUE)
all.equal(cm4[1:6], cm5[1:6])  # call and na.action won't match

u.inst <- sort(unique(tdata$inst))
temp <- matrix(0, length(u.inst), 5)
for (i in 1:length(u.inst)) {
    temp[i,] <- with(subset(tdata, inst==u.inst[i]),        
                 allpair(time, status-1, -tpred, wt))
}
aeq(temp, cm4$count)
    
eps <- 1e-6
keep <- (1:nrow(tdata))[-tfit$na.action]  # the obs that are not tossed
lmat <- matrix(0., length(keep), 5)
for (i in 1:length(keep)) {
    wt2 <- tdata$wt
    wt2[keep[i]] <- wt2[keep[i]] + eps

    test <- concordance(Surv(time, status) ~ predict(tfit) + strata(inst),
                   data=tdata, weight=wt2, group=group, reverse=TRUE,
                   keepstrata=TRUE)
    lmat[i,] <- colSums(test$count - cm4$count)/eps
}
aeq(lmat, cm4$influence, tolerance=eps)

# Check that keepstrata gives the correct sum
cm4b <- concordance(tfit, keepstrata=FALSE)
aeq(cm4b$count, colSums(cm4$count))
