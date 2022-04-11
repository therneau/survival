library(survival)
aeq <- function(x, y, ...) all.equal(as.vector(x), as.vector(y), ...)

# Check that a multi-state model, correctly set up, gives the same
# solution as a time-dependent covariate.
# This is a stronger test than mstrata: there the covariate which was mapped
#  into a state was constant, here it is time-dependent.
#
# First build the TD data set from pbcseq, with a categorical bilirubin
pbc1 <- pbcseq
pbc1$bili4 <- cut(pbc1$bili, c(0,1, 2,4, 100), 
                  c("normal", "1-2x", "2-4x", ">4"))
ptemp <- subset(pbc1, !duplicated(id))  # first row of each

pbc2 <- tmerge(ptemp[, c("id", "age", "sex")], ptemp, id,
               death= event(futime, status==2))

pbc2 <- tmerge(pbc2, pbc1, id=id, bili = tdc(day, bili),
                 bili4 = tdc(day, bili4), bstat = event(day, as.numeric(bili4)))
btemp <- with(pbc2, ifelse(death, 5, bstat))

# a row with the same starting and ending bili4 level is not an event
b2 <- ifelse(((as.numeric(pbc2$bili4)) == btemp), 0, btemp)
pbc2$bstat <- factor(b2, 0:5,
                     c("censor", "normal", "1-2x", "2-4x", ">4", "death"))
check1 <- survcheck(Surv(tstart, tstop, bstat) ~ 1, istate= bili4,
                    id = id, data=pbc2)
check1$transitions
all.equal(as.character(pbc2$bili4), as.character(check1$istate))
# the above verifies that I created the data set correctly

# Standard coxph fit with a time dependent bili4 variable.
fit1 <- coxph(Surv(tstart, tstop, death) ~ age + bili4, pbc2)

# An additive multi-state fit, where bili4 is a state
#  The three forms below should all give identical models
fit2 <- coxph(list(Surv(tstart, tstop, bstat) ~ 1,
                   c(1:4):5 ~ age / common + shared), id= id, istate=bili4,
              data=pbc2)
fit2b <- coxph(list(Surv(tstart, tstop, bstat) ~ 1,
                   1:5 + 2:5 + 3:5 + 4:5 ~ age / common + shared), 
              id= id, istate=bili4, data=pbc2)
fit2c <- coxph(list(Surv(tstart, tstop, bstat) ~ 1,
                   0:5 ~ age / common + shared), 
              id= id, istate=bili4, data=pbc2)

# Make sure the names are correct and the coefficients match
aeq(coef(fit1), coef(fit2))
aeq(names(coef(fit2)), c("age", "ph(2:5/1:5)", "ph(3:5/1:5)", "ph(4:5/1:5)"))
all.equal(coef(fit2), coef(fit2b))
all.equal(coef(fit2), coef(fit2c))

# Now a model with a separate age effect for each bilirubin group
fit3  <- coxph(Surv(tstart, tstop, death) ~ age*bili4, pbc2)
fit3b <- coxph(Surv(tstart, tstop, death) ~ bili4/age, pbc2)
fit4 <-  coxph(list(Surv(tstart, tstop, bstat) ~ 1,
                   c(1:4):5 ~ age / shared), id= id, istate=bili4,
              data=pbc2)
all.equal(fit3$loglik, fit3b$loglik)
all.equal(fit3$loglik, fit4$loglik)

# The coefficients are quite different due to different codings for dummy vars
# Unpack the interaction, first 4 coefs will be the age effect within each
#  bilirubin group
temp <- c(coef(fit3)[1] + c(0, coef(fit3)[5:7]), coef(fit3)[2:4])
names(temp)[1:4] <- c("age1", "age2", "age3", "age4")
aeq(temp, coef(fit3b)[c(4:7, 1:3)])
aeq(temp, coef(fit4))

# Third, a model with separate baseline hazards for each bili group
fit5 <- coxph(Surv(tstart, tstop, death) ~ strata(bili4)/age, pbc2,
              cluster=id)
fit6 <- coxph(list(Surv(tstart, tstop, bstat) ~ 1, 0:5 ~ age),
                   id=id, istate=bili4, pbc2)
aeq(coef(fit5), coef(fit6))
aeq(fit5$var, fit6$var)
aeq(fit5$naive.var, fit6$naive.var)


# Now for the harder part; verify that the resulting estimated survival curves
#  are correct.
# Start small with just 9 subjects, coefs fixed to values not too far from
#  the full fit.  No transtions from state 3 to death, in this subset, so
#  there is one age coefficient and 2 PH coefs
pbc3 <- subset(pbc2, id < 10)
pbc3$age <- round(pbc3$age)  # easier to do "by hand" sums
fit3 <- coxph(list(Surv(tstart, tstop, bstat) ~ 1, 
                   c(1:4):5 ~ age / common + shared), 
              id= id, istate=bili4, data=pbc3, init= c(.05, .6, 1.1), iter=0)
# a mixed p0 gives a stronger test.
surv3 <- survfit(fit3, newdata=list(age=50), p0=c(.4, .3, .2, .1, 0))

etime <- sort(unique(pbc3$tstop[pbc3$bstat != "censor"]))
# At event time 1 (182), all 9 are at risk, (3,3,2,1,0) in the initial states
atrisk <- pbc3$tstart < etime[1] & pbc3$tstop >= etime[1]  # all 9 at risk
table(pbc$bili4[atrisk])

#  One event occurs at 182, a transition from 1-2x to normal.
#  Risk scores for the non-death transitions are all exp(0) =1,
#  so the hazard matrix H will have second row of (1/3, -1/3, 0,0,0) and all
#  other rows are 0. 
with(subset(pbc3, tstop== 182), table(istate= bili4, state=bstat))

# At the next four events are from 3:4, 3:2, 2:3, and 1:2, so also have 
#  simple transtions, i.e., no covariates
#
hmat <- array(0, dim=c(5,5,6))  # first 6 hazard matrices, start with 3,3,2,1,0
hmat[2,1,1] <- 1/3; hmat[2,2,1] <- -1/3   # new count= 4,2,2,1,0
hmat[3,4,2] <- 1/2; hmat[3,3,2] <- -1/2   # new count= 4,2,1,2,0
hmat[3,2,3] <- 1  ; hmat[3,3,3] <- -1     # new count= 4,3,0,2,0
hmat[2,3,4] <- 1/3; hmat[2,2,4] <- -1/3   # new count= 4,2,1,2,0
hmat[1,2,5] <- 1/4; hmat[1,1,5] <- -1/4   # new count= 3,3,1,2,0

# Event 6 is a transition from state 4 to death
# For the shared hazard, the denominator is all those in states 1,2, or 4.
#
atrisk <- subset(pbc3, tstart < etime[6] & tstop >= etime[6] & bili4 != '2-4x')
eta <- with(atrisk, .05*(age-50) + .6*(bili4=="1-2x") + 1.1*(bili4 == ">4x"))
basehaz <- 1/sum(exp(eta))
hmat[1,5,6] <- basehaz;            hmat[1,1,6] <- -basehaz
hmat[2,5,6] <- basehaz * exp(.7);  hmat[2,2,6] <- -basehaz*exp(.6)
hmat[4,5,6] <- basehaz * exp(1.1); hmat[4,4,6] <- -basehaz*exp(1.1)

tmat <- hmat  # transition matrices
pstate <- matrix((4:0)/10, nrow=1)
for (i in 1:6) {
    tmat[,,i] <- as.matrix(expm(hmat[,,i]))
    pstate <- rbind(pstate, pstate[i,]%*% tmat[,,i])
}

test3 <- with(p3, mysurv(tstart, tstop, bstat, bili4, p0= (4:0)/10, fit3,
                         50))

#
# Creating risk sets one by one like this if far too slow for the production
#  code; which does a lot more preemptive bookkeeping.
#
# The linear predictors and the strata map are pulled from the fit.
#  
mysurv <- function(t1, t2, state, istate, p0, fit, baseline, debug=0) {
    if (!inherits(fit, 'coxphms')) stop("invalid fit")
    smap <- fit$stratum_map
    from <- as.numeric(sub(":.*$", "", colnames(smap)))
    to   <- as.numeric(sub("^.*:", "", colnames(smap)))
    shared <- duplicated(smap[1,])
    nshare <- sum(shared)
    bcoef <- rep(1, ncol(smap))   # coefficients for shared baseline
    beta <- coef(fit, matrix=TRUE)
    if (nshare >0) {
        # coefficients for shared baseline will be the last nshare of them
        i <- seq(length=nshare, to=length(fit$coefficients))
        bcoef[shared] <- exp(fit$coefficients[i])
        # remove shared coef rows from beta
        phrow <- apply(fit$cmap, 1, function(x) any(x %in% i))
        beta <- beta[!phrow,]
    }
    baserisk <- baseline %*% beta

    # Make the values for istate and state match the 1:2, etc of the fit,
    #  i.e., the order of fit$states
    istate <- factor(as.character(istate), fit$states)
    state  <- factor(as.character(state),  fit$states) #censor= NA
    nstate <- length(fit$states)

    # set up output
    ntran <- ncol(smap)    # number of transitions
    utime <- sort(unique(t2[!is.na(state)])) # unique event times
    ntime <- length(utime)
    tmat <- matrix(0, nstate, nstate)  # transtion matrix at this time point
    pmat <- diag(nstate)               # product of transitions
    nrisk <- matrix(0., ntime, nstate) #number at risk
    wtrisk<- matrix(0., ntime, ntran)  # weighted number per transtion 
    nevent <- matrix(0L, ntime, nstate)  # number of events of each type
    pstate <- matrix(0L, ntime, nstate)  # probability in state
    hmat <- matrix(0., nstate, nstate)   # working matrix of hazards
    rwt <- exp(fit$linear.predictor - rep(baserisk, each=length(t1)))
    if (nrow(rwt) != length(t1)) stop("mismatched time vector")

    for (i in 1:ntime) {
        atrisk <- (t1 < utime[i] &  utime[i] <= t2) # risk set at this time
        event <- which(utime[i] == t2)  # potential events, end at this time
        nrisk[i,]  <- c(table(istate[atrisk]))   # number at risk in each state
        nevent[i,] <- c(table(state[event]))
        # The linear predictor and hence the number at risk is different for
        #  every transition. Also, so will not be at risk for the transition.
        #
        for (k in 1:ntran) { 
            atrisk2 <- (atrisk & (as.numeric(istate) == from[k]))
            wtrisk[i,k] <- sum(rwt[atrisk2,k])
            }
        dtemp <- table(istate[event], state[event]) #censors don't count

        # fill in hmat, one hazard at a time
        hmat <- 0*hmat 
        for (j in unique(smap)) {
            # for each baseline hazard
            k <- which(smap == j)  # transitons that share this hazard
            deaths <- sum(dtemp[cbind(from[k], to[k])])     # total events
            if (deaths==0) hmat[cbind(from[k], to[k])] <- 0 # avoid 0/0
            else {
                hazard <- deaths/ sum(wtrisk[i, k] * bcoef[k]) #shared baseline
                hmat[cbind(from[k], to[k])] <- hazard * bcoef[k] # PH
            }
        }
        diag(hmat) <- diag(hmat) - rowSums(hmat)   # rows sum to zero
        tmat <- as.matrix(expm(hmat))          # transtion matrix
        if (debug>0) browser()
        pmat <- pmat %*% tmat
        pstate[i,] <- drop(p0 %*% pmat)
    }
    list(time=utime, nrisk=nrisk, nevent=nevent, pstate=pstate,  
         wtrisk= wtrisk, P=pmat)
}

dummy <- data.frame(age = 45)
surv2 <- survfit(fit2, newdata=dummy, p0=c(1,0,0,0,0))

# risk score for each observation, relative to a 45 year old
rwt <- exp(predict(fit2) - 45*fit2$coef[1])
test2 <- with(pbc2, mysurv(tstart, tstop, bstat, bili4, rwt= rwt, 
                  p0= c(1,0,0,0,0), smap= fit2$stratum_map, fit2$coef[2:4]))
if (FALSE){
    # for testing
    plot(surv2, col=1:5, lwd=2)
    matpoints(test2$time, test2$pstate, col=1:5, pch='o')
}

# with risk scores of all 1
ones <- rep(1, nrow(pbc2))
rtest <-  with(pbc2, mysurv(tstart, tstop, bstat, bili4, ones, 
                            p0=c(1,0,0,0,0), test2$stratum_map))
test3 <- survfit(Surv(tstart, tstop, bstat) ~1, pbc2, id=id, istate=bili4,
                 p0 = c(1,0,0,0,0))


# smaller, only 2 states transtion to death
p3 <- subset(pbc2, id< 10)
fit3 <- coxph(list(Surv(tstart, tstop, bstat) ~ 1,
                   c(1:4):5 ~ age / common + shared), id= id, istate=bili4,
              data=p3, iter=1)
surv3 <- survfit(fit3, newdata= dummy, p0=c(1,0,0,0,0))

rwt <- exp(predict(fit3) - 45*fit3$coef[1])
test3 <- with(p3, mysurv(tstart, tstop, bstat, bili4, p0= c(1,0,0,0,0), fit3,
                         45))
