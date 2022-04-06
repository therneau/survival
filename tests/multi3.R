library(survival)
aeq <- function(x, y, ...) all.equal(as.vector(x), as.vector(y), ...)

# Check that a multi-state model, correctly set up, gives the same
# solution as a time-dependent covariate.
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

# Time dependent fit
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

# unpack the interaction, first 4 coefs will be the age effect within each
#  bilirubin group
temp <- c(coef(fit3)[1] + c(0, coef(fit3)[5:7]), coef(fit3)[2:4])
names(temp)[1:4] <- c("age1", "age2", "age3", "age4")
aeq(temp, coef(fit3b)[c(4:7, 1:3)])
aeq(temp, coef(fit4))

# Third, a model with separate baseline hazards
fit5 <- coxph(Surv(tstart, tstop, death) ~ strata(bili4)/age, pbc2)
fit6 <- coxph(list(Surv(tstart, tstop, bstat) ~ 1, 0:5 ~ age),
                   id=id, istate=bili4, pbc2)
aeq(coef(fit5), coef(fit6))

# Now for the hard part; verify that the resulting estimated survival curves
#  are correct
# First fit1/fit2

# Compute the curve by hand
# (Creating risk sets one by one like this if far too slow for the production
#  code; which does a lot more preemptive bookkeeping.)
#
# The last arguements are the risk score eta for each observation, which is
#  constant for any subject, and the strata map, and shared strata coefs
coxsurve <- function(t1, t2, state, istate, rwt, p0, smap, scoef) {
    if (missing(smap)) stop("a strata map is required")
    # set up shared baseline coefficients = bcoef
    temp <- colnames(smap)
    from <- as.numeric(sub(":.*$", "", temp))
    to   <- as.numeric(sub("^.*:", "", temp))
    smap <- smap[1,]  # we only need the baseline part
    shared <- duplicated(smap)
    nshare <- sum(shared)
    bcoef <- rep(1, length(smap))   # coefficients for shared baseline
    if (nshare >0) {
        if (missing(scoef) || length(scoef) != nshare) stop("invalid scoef")
        bcoef[shared] <- exp(scoef)
    }

    # set up istate and state to have matching levels, state also has 'censored'
    cstate <- as.numeric(istate)
    allstate <- unique(c(levels(state), levels(istate)))
    nstate <- length(allstate) -1
    # make istate have all the levels
    if (length(levels(istate)) != nstate) 
        istate <- factor(as.numeric(istate), 1:nstate, allstate[-1])
    # and state too
    if (length(levels(state)) != nstate +1)
        state <- factor(as.numeric(state), 1:(nstate+1), allstate)

    # set up output
    utime <- sort(unique(t2[state != levels(state)[1]]))
    ntime <- length(utime)
    tmat <- matrix(0, nstate, nstate) # transtion matrix at this time point
    pmat <- diag(nstate)              # product of transitions
    nrisk <- wtrisk <- matrix(0., ntime, nstate) #number at risk, weighted nrisk
    nevent <- matrix(0L, ntime, nstate)  # number of events of each type
    pstate <- matrix(0L, ntime, nstate)  # probability in state
    hmat <- matrix(0., nstate, nstate)   # working matrix of hazards
    for (i in 1:ntime) {
        atrisk <- which(t1 < utime[i] &  utime[i] <= t2) # risk set at this time
        event <- which(utime[i] == t2)  # events at this time
        nrisk[i,]  <- c(table(istate[atrisk]))   # number at risk in each state
        temp  <- c(tapply(rwt[atrisk], istate[atrisk], sum))
        wtrisk[i,] <- ifelse(is.na(temp), 0, temp)  # tapply creates NA
        dtemp <- table(istate[event], state[event])[,-1]
        for (j in unique(smap)) {
            k <- which(smap == j)
            deaths <- sum(dtemp[cbind(from[k], to[k])]) 
            if (deaths==0) hmat[cbind(from[k], to[k])] <- 0
            else  hmat[cbind(from[k], to[k])] <- deaths/
                      sum(wtrisk[i, from[k]] * bcoef[k])
        }
        diag(hmat) <- diag(hmat) - rowSums(hmat)   # rows sum to zero
        pmat <- pmat %*% expm(hmat)
        pstate[i,] <- p0 <- drop(p0 %*% expm(hmat))
    }
    list(time=utime, nrisk=nrisk, nevent=nevent, pstate=pstate,  
         wtrisk= wtrisk, P=pmat)
}

dummy <- data.frame(age = 45)
surv2 <- survfit(fit2, newdata=dummy, p0=c(1,0,0,0,0))

# risk score for each subject, at event times, relative to a 45 year old
rwt <- exp(predict(fit2) - 45*fit2$coef[1])
test2 <- with(pbc2, coxsurve(tstart, tstop, bstat, bili4, rwt= rwt, 
                  p0= c(1,0,0,0,0), smap= fit2$stratum_map, fit2$coef[2:4]))

ones <- rep(1, nrow(pbc2))
rtest <-  with(pbc2, handcurve(tstart, tstop, bstat, bili4, ones, c(1,0,0,0,0)))
test3 <- survfit(Surv(tstart, tstop, bstat) ~1, pbc2, id=id, istate=bili4,
                 p0 = c(1,0,0,0,0))
plot(test3, col=1:5)
matpoints(rtest$time, rtest$pstate, pch=1, col=1:5)
