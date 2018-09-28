library(survival)
options(na.action=na.exclude) # preserve missings
options(contrasts=c('contr.treatment', 'contr.poly')) #ensure constrast type

#
# Tests for the condordance function. 
#
aeq <- function(x,y, ...) all.equal(as.vector(x), as.vector(y), ...)

grank <- function(x, time, grp, wt) {
    if (missing(wt) || length(wt)==0) 2*unlist(tapply(x, grp, rank))
    else {
        z <- double(length(x))
        for (i in unique(grp)) {
            indx <- which(grp==i)
            temp <- tapply(wt[indx], x[indx], sum)
            temp <- temp/2  + c(0, cumsum(temp)[-length(temp)])
            z[indx] <- temp[match(x[indx], names(temp))]
        }
        z
    }
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

# PH variance by brute force
phvar <- function(time, status, x, wt) {
    if (missing(wt)) wt <- rep(1, length(x))
    zmat <- wt* outer(x, x, function(x, y) sign(x-y))
    z2 <- sapply(which(status==1), function(i) {
            atrisk <- (time >= time[i])
            zscore <- colSums(zmat[atrisk,, drop=FALSE])
            c(zscore[i], sum((wt*zscore^2)[atrisk])/sum(wt[atrisk]))
    })
    rowSums(z2 * rep(wt[status==1], each=2))  # Cox score stat, var of score
}
                      
tdata <- aml[aml$x=='Maintained', c("time", "status")]
tdata$x <- c(1,6,2,7,3,7,3,8,4,4,5)
tdata$wt <- c(1,2,3,2,1,2,3,4,3,2,1)

fit <- concordance(Surv(time, status) ~x, tdata, influence=2)
aeq(fit$count, with(tdata, allpair(time, status, x)))
aeq(fit$influence, with(tdata, leverage(time, status, x)))
npair <- sum(fit$count[1:3])
aeq(c(fit$count[1]-fit$count[2], 4*npair^2*fit$var[2]), 
    with(tdata, phvar(time, status, x)))

# verify that the phvar function is correct by fitting a time-dependent
#  Cox model
cfit <- coxph(Surv(time, status) ~ tt(x), tdata, tt=grank, method='breslow',
              iter=0, x=T, weights=wt)
cdt <- coxph.detail(cfit)
aeq(c(-sum(cdt$score), sum(cdt$imat)),  with(tdata, phvar(time, status, x, wt)))

# Weighted
fit <- concordance(Surv(time, status) ~x, tdata, influence=2, weights=wt)
aeq(fit$count, with(tdata, allpair(time, status, x, wt)))
aeq(fit$influence, with(tdata, leverage(time, status, x, wt)))
npair <- sum(fit$count[1:3])
aeq(c(fit$count[1]-fit$count[2], 4*npair^2*fit$var[2]), 
    with(tdata, phvar(time, status, x, wt)))


# Test 2: Lots of ties
tempy <- Surv(c(1,2,2,2,3,4,4,4,5,2), c(1,0,1,0,1,0,1,1,0,1))
tempx <- c(5,5,4,4,3,3,7,6,5,4)
fit2 <- concordance(tempy ~ tempx, influence=2)
aeq(fit2$count, allpair(tempy[,1], tempy[,2], tempx))
aeq(fit2$influence, leverage(tempy[,1], tempy[,2], tempx))
npair <- sum(fit2$count[1:3])
aeq(4*npair^2*fit2$var[2], phvar(tempy[,1], tempy[,2], tempx)[2])

# Bigger data
cox3 <- coxph(Surv(time, status) ~ age + sex + ph.ecog, lung)
fit3 <- concordance(Surv(time, status) ~ predict(cox3), lung, influence=2)
tdata <- na.omit(lung[,c('time', 'status', 'age', 'sex', 'ph.ecog')])
aeq(fit3$count, allpair(lung$time, lung$status-1,predict(cox3)))
aeq(fit3$influence, leverage(lung$time, lung$status-1,predict(cox3)))


# More ties
fit4 <- concordance(Surv(time, status) ~ ph.ecog, lung, reverse=TRUE)
aeq(fit4$count, allpair(lung$ph.ecog, lung$time, lung$status-1))
aeq(fit4$count, c(8392, 4258, 7137, 28))
cfit4 <- coxph(Surv(time, status) ~ tt(ph.ecog), lung, 
               iter=0, method='breslow', tt=grank)
aeq(4/cfit4$var, fit4$stats[5]^2)

# Case weights
fit5 <- concordance(Surv(time, status) ~ x, tdata, weight=wt)
fit6 <- concordance(Surv(time, status) ~x, tdata[rep(1:11,tdata$wt),])
aeq(fit5$count, with(tdata, allpair(x, time, status, wt)))
aeq(fit5$count, c(91, 70, 7, 0))  # checked by hand
aeq(fit5$count[1:3], fit6$count[1:3])  #spurious "tied on time" values, ignore
aeq(fit5$std, fit6$std)
cfit5 <- coxph(Surv(time, status) ~ tt(y), tdata, weight=wt, 
               iter=0, method='breslow', tt=grank2)
cfit6 <- coxph(Surv(time, status) ~ tt(y), tdata[rep(1:11,tdata$wt),], 
               iter=0, method='breslow', tt=grank)
aeq(4/cfit6$var, fit6$stats[5]^2)
aeq(cfit5$var, cfit6$var)

# Start, stop simplest cases
fit7 <- concordance(Surv(rep(0,11), time, status) ~ y, tdata)
aeq(fit7$stats, fit$stats)
aeq(fit7$std.err, fit$std.err)
fit7 <- concordance(Surv(rep(0,11), time, status) ~ y, tdata, weight=wt)
aeq(fit5$stats, fit7$stats)

# Multiple intervals for some, but same risk sets as tdata
tdata2 <- data.frame(time1=c(0,3, 5,  6,7,   0,  4,17,  7,  0,16,  2,  0, 
                             0,9, 5),
                     time2=c(3,9, 13, 7,13, 18, 17,23, 28, 16,31, 34, 45, 
                             9,48, 60),
                     status=c(0,1, 1, 0,0,  1,  0,1, 0, 0,1, 1, 0, 0,1, 0),
                     x = c(1,1, 6, 2,2, 7, 3,3, 7, 3,3, 8, 4, 4,4, 5),
                     wt= c(1,1, 2, 3,3, 2, 1,1, 2, 3,3, 4, 3, 2,2, 1))
fit8 <- concordance(Surv(time1, time2, status) ~x, tdata2, weight=wt)
aeq(fit5$stats, fit8$stats)
aeq(fit5$std.err, fit8$std.err)
cfit8 <- coxph(Surv(time1, time2, status) ~ tt(x), tdata2, weight=wt, 
               iter=0, method='breslow', tt=grank2)
aeq(4/cfit8$var, fit8$stats[5]^2)
aeq(fit8$stats[5]/(2*sum(fit8$stats[1:3])), fit8$std.err)

# Stratified
tdata3 <- data.frame(time1=c(tdata2$time1, rep(0, nrow(lung))),
                     time2=c(tdata2$time2, lung$time),
                     status = c(tdata2$status, lung$status -1),
                     x = c(tdata2$x, lung$ph.ecog),
                     wt= c(tdata2$wt, rep(1, nrow(lung))),
                     grp=rep(1:2, c(nrow(tdata2), nrow(lung))))
fit9 <- concordance(Surv(time1, time2, status) ~x + strata(grp),
                        data=tdata3, weight=wt)
aeq(fit9$stats[1,], fit5$stats)
aeq(fit9$stats[2,], fit4$stats)
