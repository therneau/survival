#
# Check out the survfit routine on the simple AML data set.
#  The leverage validation makes use of the fact that when all
#  weights are 1 and there is 1 obs per subject, the IJ variance is
#  equal to the Greenwood.
# There are 8 choices in the C code:  Nelson-Aalen or Fleming-Harrington
#  estimate of cumulative hazard,  KM or exp(cumhaz) estimate of survival,
#  regular or robust variance.  This tries to exercise them all.

library(survival)
aeq <- function(x, y, ...) all.equal(as.vector(x), as.vector(y), ...)

set.seed(1953)  # used only to reorder the data
adata <- aml
adata$id <- sample(LETTERS, nrow(aml)) # labels are not in time or data order
adata <- adata[sample(1:nrow(aml), nrow(aml)),] # data is unordered
adata$wt <- sample((2:30)/10, nrow(aml))  # non-integer weights

byhand <- function(time, status, weights, id) {
    # for a single curve
    utime <- sort(unique(time))
    ntime <- length(utime)
    n <- length(time)
    if (missing(weights)) weights <- rep(1.0, n)
    if (missing(id)) id <- seq(along=time)

    uid <- sort(unique(id))
    nid <- length(uid)
    id <- match(id, uid)  # change it to 1:nid

    n.risk <- n.event <- surv <- cumhaz <- double(ntime)
    KM <- 1; nelson <-0; 
    kvar <- 0; hvar<-0;
        
    U <- matrix(0, nid, 2)  # the two robust influence estimates
    V <- matrix(0, ntime, 4)  # variances
    usave <- array(0., dim=c(nid, 2, ntime))
    estimate <- matrix(0, ntime, 2)

    for (i in 1:ntime) {
        atrisk <- (time >= utime[i])
        n.risk[i] <- sum(weights[atrisk])
        deaths <- (time==utime[i] & status==1)                 
        n.event[i] <- sum(weights[deaths])
 
        haz <- n.event[i]/n.risk[i]
        dhaz <- (ifelse(deaths,1,0) - ifelse(atrisk, haz, 0))/n.risk[i]
        U[,1] <- U[,1]*(1-haz) - KM*tapply(dhaz*weights, id, sum)
        V[i,1] <- sum(U[,1]^2)
            
        U[,2] <- U[,2] + tapply(dhaz* weights, id, sum) #result in 'id' order
        V[i,2] <- sum(U[,2]^2)
        usave[,,i] <- U

        if (n.event[i] >0 ) {
            KM <- KM*(1-haz)
            nelson <- nelson + haz
            kvar <- kvar + n.event[i]/(n.risk[i] * (n.risk[i] - n.event[i]))
            hvar <- hvar + n.event[i]/(n.risk[i]^2)
        }
            
        V[i,3] <- kvar   # var of log(S)
        V[i,4] <- hvar
        estimate[i,] <- c(KM, nelson)
        }
    dimnames(usave) <- list(uid, c("KM", "chaz"), utime)
    list(time=utime, n.risk=n.risk, n.event=n.event, estimate=estimate,
         std = sqrt(V), influence=usave)
}

# the byhand function can only handle one group at a time
true1a <- with(subset(adata, x=="Maintained"), byhand(time, status, id=id))
true1b <- with(subset(adata, x!="Maintained"), byhand(time, status, id=id))

# The Greenwood and IJ estimates agree, except for a last point with
#  variance of zero.  These next few lines verify the byhand() function
aeq(true1a$std[,1], true1a$estimate[,1]*true1a$std[,3])
aeq(true1b$std[1:9,1], true1b$estimate[1:9,1]*true1b$std[1:9,3])
aeq(true1b$std[10,1], 0)   # variance of zero for jackknife
!is.finite(true1b$std[10,3])   # Inf for Greenwood
temp <- with(subset(adata, x=="Maintained"), byhand(time, status, id=id,
                                                    weights=rep(3,11)))
aeq(temp$std[,1:2], true1a$std[,1:2])  # IJ estimates should be invariant

# fit1 uses the standard formulas: NA hazard, KM survival
fit1 <- survfit(Surv(time, status) ~ x, data=adata)
aeq(fit1$surv, c(true1a$estimate[,1], true1b$estimate[,1]))
aeq(fit1$cumhaz, c(true1a$estimate[,2], true1b$estimate[,2]))
aeq(fit1$std.err, c(true1a$std[,3], true1b$std[,3]))
aeq(fit1$std.chaz, c(true1a$std[,4], true1b$std[,4]))
aeq(fit1$n.risk, c(true1a$n.risk, true1b$n.risk))
aeq(fit1$n.event, c(true1a$n.event, true1b$n.event))
fit1$logse   # logse should be TRUE

# fit2 will use the IJ method
fit2 <- survfit(Surv(time, status) ~ x, data=adata, id=id, robust=TRUE)
aeq(fit2$surv, c(true1a$estimate[,1], true1b$estimate[,1]))
aeq(fit2$cumhaz, c(true1a$estimate[,2], true1b$estimate[,2]))
aeq(fit2$std.err, c(true1a$std[,1], true1b$std[,1]))
aeq(fit2$std.chaz, c(true1a$std[,2], true1b$std[,2]))
aeq(fit2$n.risk, c(true1a$n.risk, true1b$n.risk))
aeq(fit2$n.event, c(true1a$n.event, true1b$n.event))
!fit2$logse  # logse should be FALSE

# look at the leverage values
fit3 <- survfit(Surv(time, status) ~ x, data=adata, id=id, influence=3)
aeq(fit3$influence.surv[[1]], true1a$influence[,1,])
aeq(fit3$influence.surv[[2]], true1b$influence[,1,])
aeq(fit3$influence.chaz[[1]], true1a$influence[,2,])
aeq(fit3$influence.chaz[[2]], true1b$influence[,2,])

# compute the influence by brute force
tdata <- subset(adata, x != "Maintained")
tdata <- tdata[order(tdata$id),]   # easier to compare if it's in order
eps <- 1e-8
imat1 <- imat2 <-  matrix(0., 12, 10)
t1 <- survfit(Surv(time, status) ~x, data=tdata) 
for (i in 1:12) {
    wtemp <- rep(1.0, 12)
    wtemp[i] <- 1 + eps
    tfit <-survfit(Surv(time, status) ~x, data=tdata, weights=wtemp) 
    imat2[i,] <- (tfit$cumhaz - t1$cumhaz)/eps
    imat1[i,] <- (tfit$surv - t1$surv)/eps
}
aeq(imat1, true1b$influence[,1,], tol= sqrt(eps))
aeq(imat2, true1b$influence[,2,], tol= sqrt(eps))

# Repeat using the Nelson-Aalen hazard and exp(NA) for survival
fit1 <- survfit(Surv(time, status) ~ x, adata, stype=2)
aeq(fit1$surv, exp(-c(true1a$estimate[,2], true1b$estimate[,2])))
aeq(fit1$cumhaz, c(true1a$estimate[,2], true1b$estimate[,2]))
aeq(fit1$std.err, c(true1a$std[,4], true1b$std[,4]))
aeq(fit1$std.chaz, c(true1a$std[,4], true1b$std[,4]))
aeq(fit1$n.risk, c(true1a$n.risk, true1b$n.risk))

# Nelson-Aalen + exp() surv, along with IJ variance
fit2 <- survfit(Surv(time, status) ~ x, data=adata, id=id, stype=2,
                influence=3)
aeq(fit2$surv, exp(-c(true1a$estimate[,2], true1b$estimate[,2])))
aeq(fit2$cumhaz, c(true1a$estimate[,2], true1b$estimate[,2]))
aeq(fit2$std.err, c(true1a$std[,2], true1b$std[,2]))
aeq(fit2$std.chaz, c(true1a$std[,2], true1b$std[,2]))
aeq(fit2$n.risk, c(true1a$n.risk, true1b$n.risk))
aeq(fit2$influence.chaz[[1]], true1a$influence[,2,])
aeq(fit2$influence.chaz[[2]], true1b$influence[,2,])
aeq(fit2$influence.surv[[2]], -true1b$influence[,2,]%*% diag(fit2[2]$surv))

# Cumulative hazard is the same for fit1 and fit2
all.equal(fit2$influence.chaz, fit3$influence.chaz)

# Weighted fits
true2a <- with(subset(adata, x=="Maintained"), byhand(time, status, id=id,
                                                      weights= wt))
true2b <- with(subset(adata, x!="Maintained"), byhand(time, status, id=id,
                                                      weights=wt))
fit3 <- survfit(Surv(time, status) ~ x, data=adata, id=id, weights=wt,
                 influence=TRUE)  
aeq(fit3$influence.surv[[1]], true2a$influence[,1,])
aeq(fit3$influence.surv[[2]], true2b$influence[,1,])
aeq(fit3$influence.chaz[[1]], true2a$influence[,2,])
aeq(fit3$influence.chaz[[2]], true2b$influence[,2,])
aeq(fit3$surv, c(true2a$estimate[,1], true2b$estimate[,1]))
aeq(fit3$cumhaz, c(true2a$estimate[,2], true2b$estimate[,2]))
aeq(fit3$std.err, c(true2a$std[,1], true2b$std[,1]))
aeq(fit3$std.chaz, c(true2a$std[,2], true2b$std[,2]))
aeq(fit3$n.risk, c(true2a$n.risk, true2b$n.risk))
aeq(fit3$n.event, c(true2a$n.event, true2b$n.event))

# Different survival, same hazard
fit3b <- survfit(Surv(time, status) ~ x, data=adata, id=id, weights=wt,
                 influence=2, stype=2) 
temp <- c("n", "time", "cumhaz", "std.chaz", "influence.chaz", "n.risk",
          "n.event")
aeq(unclass(fit3b)[temp], unclass(fit3)[temp])  # unclass avoids [.survfit
aeq(fit3b$surv, exp(-c(true2a$estimate[,2], true2b$estimate[,2])))
aeq(fit3b$std.err, fit3b$std.chaz)
aeq(fit3b$logse, FALSE)
aeq(fit3b$n.risk, c(true2a$n.risk, true2b$n.risk))
aeq(fit3b$n.event, c(true2a$n.event, true2b$n.event))

# The grouped jackknife
group <- rep("", nrow(adata))
temp <- table(adata$x)
group[adata$x == "Maintained"] <- rep(letters[1:4], length=temp[1])
group[adata$x != "Maintained"] <- rep(letters[4:7], length=temp[2])
adata$group <- group
fit4 <-  survfit(Surv(time, status) ~ x, data=adata, id=id, weights=wt,
                 influence=TRUE, cluster=group)
g1 <- adata$group[match(rownames(true2a$influence[,1,]), adata$id)]
g2 <- adata$group[match(rownames(true2b$influence[,1,]), adata$id)] 
aeq(fit4$influence.surv[[1]], rowsum(true2a$influence[,1,], g1))
aeq(fit4$influence.surv[[2]], rowsum(true2b$influence[,1,], g2))
aeq(fit4$influence.chaz[[1]], rowsum(true2a$influence[,2,], g1))
aeq(fit4$influence.chaz[[2]], rowsum(true2b$influence[,2,], g2))

aeq(c(colSums(fit4$influence.surv[[1]]^2), colSums(fit4$influence.surv[[2]]^2)),
    fit4$std.err^2)
aeq(c(colSums(fit4$influence.chaz[[1]]^2), colSums(fit4$influence.chaz[[2]]^2)),
    fit4$std.chaz^2)

# The Fleming-Harrington is a more complex formula.  Start with weights of
#   1.
fit5 <- survfit(Surv(time, status) ~x, adata, ctype=2)
nrisk <- c(11,10,8,7, 5,4,2, 12, 11, 10, 9, 8, 6:1)
chaz <- c(cumsum(1/nrisk[1:7])[c(1:4,4, 5,6,6,7,7)], 
          cumsum(1/nrisk[8:18])[c(2,4,5,5,6:11)])
aeq(fit5$cumhaz, chaz)
aeq(fit5$std.chaz, sqrt(c(cumsum(1/nrisk[1:7]^2)[c(1:4,4, 5,6,6,7,7)], 
                          cumsum(1/nrisk[8:18]^2)[c(2,4,5,5,6:11)])))

# We can compute the FH using a fake data set where each tie is spread out
#  over a set of fake times.
# 
fh <- function(time, status, weights, id) {
    counts <- table(time, status)
    utime <-  sort(unique(time))
    tied <- counts[,2] > 1

    if (missing(weights)) weights <- rep(1.0, length(time))
    if (missing(id))  id <- 1:length(time)

    # build the expanded data set
    delta <- min(diff(utime))/(2*max(counts[,2]))
    efun <- function(x) {
        who <- which(time==x & status==1)
        ntie <- length(who)
        data.frame(time = rep(x - (1:ntie -1)*delta, each=ntie),
                   id = rep(id[who], ntie),
                   status = rep(1, ntie^2),
                   weight = rep(weights[who]/ntie, ntie),
                   stringsAsFactors=FALSE
                   )
    }

    temp <- do.call(rbind, lapply(utime[tied], efun))
    notie <- (status==0 | !(time %in% utime[tied]))

    bfit <- byhand(time = c(time[notie], temp$time), 
                   status = c(status[notie], temp$status),
                   id = c(id[notie], temp$id),
                   weights = c(weights[notie], temp$weight)
                   )
    keep <- match(utime, bfit$time)  # the real time points

    list(time=bfit$time[keep], 
         n.risk=bfit$n.risk[keep - pmax(0, counts[,2]-1)],
         n.event = bfit$n.event[keep]* counts[,2],  
         estimate=bfit$estimate[keep,],
         std = bfit$std[keep,], influence=bfit$influence[,,keep])
}

# Case weights
true6a <- with(subset(adata, x=="Maintained"), fh(time, status, wt, id))
true6b <- with(subset(adata, x!="Maintained"), fh(time, status, wt, id))

fit6 <- survfit(Surv(time, status) ~ x, weight=wt, data=adata, stype=2, 
                ctype=2, robust=FALSE)
aeq(fit6$cumhaz, c(true6a$estimate[,2], true6b$estimate[,2]))
aeq(fit6$surv, exp(-c(true6a$estimate[,2], true6b$estimate[,2])))
aeq(fit6$std.chaz, c(true6a$std[,4], true6b$std[,4]))
aeq(fit6$n.risk, c(true6a$n.risk, true6b$n.risk))
aeq(fit6$n.event, c(true6a$n.event, true6b$n.event))

# Robust variance
fit7 <- survfit(Surv(time, status) ~ x, weight=wt, data=adata, stype=2,ctype=2, 
                id=id, influence=2, robust=TRUE)
aeq(fit7$cumhaz, c(true6a$estimate[,2], true6b$estimate[,2]))
aeq(fit7$surv, exp(-c(true6a$estimate[,2], true6b$estimate[,2])))
aeq(fit7$std.chaz, c(true6a$std[,2], true6b$std[,2]))
aeq(fit7$n.risk, c(true6a$n.risk, true6b$n.risk))
aeq(fit7$n.event, c(true6a$n.event, true6b$n.event))
aeq(fit7$influence.chaz[[1]], true6a$influence[,2,])
aeq(fit7$influence.chaz[[2]], true6b$influence[,2,])

# compute the influence by brute force
tdata <- subset(adata, x != "Maintained")
tdata <- tdata[order(tdata$id),]
eps <- 1e-8
imat <- matrix(0., 12, 10)
t1 <- survfit(Surv(time, status) ~x, data=tdata, ctype=2, weights=wt) 
for (i in 1:12) {
    wtemp <- tdata$wt
    wtemp[i] <- wtemp[i] + eps
    tfit <-survfit(Surv(time, status) ~x, data=tdata, ctype=2, 
              weights=wtemp)
    imat[i,] <- tdata$wt[i] * (tfit$cumhaz - t1$cumhaz)/eps
}
aeq(fit7$influence.chaz[[2]], imat, tol=sqrt(eps))

#
# verify that the times and scale arguments work as expected.  They
#  are in the summary and print.survfit functions.
#
s1 <- summary(fit1, scale=1)
s2 <- summary(fit1, scale=2)
aeq(s1$time/2, s2$time)  #times change
aeq(s1$surv, s2$surv)
tscale <- rep(c(1,1,1,1, 2,2,2,2,2), each=2)  
aeq(s1$table, s2$table *tscale)

s3 <- summary(fit1, scale=1, times=c(9, 18, 23, 33, 34))
s4 <- summary(fit1, scale=2, times=c(9, 18, 23, 33, 34))
aeq(s3$time, s4$time*2)
aeq(s3$surv, s4$surv)

print(fit1, rmean='common')
print(fit1, rmean='common', scale=2)
