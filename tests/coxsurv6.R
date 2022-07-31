library(survival)
aeq <- function(x, y, ...) all.equal(as.vector(x), as.vector(y), ...)

# Test the survival curve for a fit with shared hazards.
# Use the pbcseq data set, and turn bilirubin into a time-dependent state with
#  4 levels, and a shared baseline hazard for the 4 transitions to death.
# The subtlety is that coefficients for a shared (proportional) baseline hazard
#  are attached to a state, not to an observation. 
# (A bilirubin value of <1 is normal.)
pbc1 <- pbcseq
pbc1$bili4 <- cut(pbc1$bili, c(0,1, 2,4, 100), 
                  c("normal", "1-2", "2-4", ">4"))
ptemp <- subset(pbc1, !duplicated(id))  # first row of each

pbc2 <- tmerge(ptemp[, c("id", "age", "sex")], ptemp, id,
               death= event(futime, status==2))

pbc2 <- tmerge(pbc2, pbc1, id=id, bili = tdc(day, bili),
                 bili4 = tdc(day, bili4), bstat = event(day, as.numeric(bili4)))
btemp <- with(pbc2, ifelse(death, 5, bstat))

# a row with the same starting and ending bili4 level is not an event
b2 <- ifelse(((as.numeric(pbc2$bili4)) == btemp), 0, btemp)
pbc2$bstat <- factor(b2, 0:5,
                     c("censor", "normal", "1-2", "2-4", ">4", "death"))
check1 <- survcheck(Surv(tstart, tstop, bstat) ~ 1, istate= bili4,
                    id = id, data=pbc2)
check1$transitions

fit2 <- coxph(list(Surv(tstart, tstop, bstat) ~ 1,
                   c(1:4):5 ~ age / common + shared), id= id, istate=bili4,
              data=pbc2)

# Before we tackle fit2, start small with just 9 subjects, coefs fixed to 
# simple values to make hand computation easier.  There are no transitions
# from state 3 to death in this subset, so there is one age coefficient and 
# 2 PH coefs.  
pbc3 <- subset(pbc2, id < 10)
pbc3$age <- round(pbc3$age)  # easier to do "by hand" sums
fit3 <- coxph(list(Surv(tstart, tstop, bstat) ~ 1, 
                   c(1:4):5 ~ age / common + shared),  x=TRUE,
              id= id, istate=bili4, data=pbc3, init= c(.05, .6, 1.1), iter=0)
# a mixed p0 gives a stronger test than our usual (1, 0,0,0,0)
surv3 <- survfit(fit3, newdata=list(age=50), p0=c(.4, .3, .2, .1, 0))

etime <- sort(unique(pbc3$tstop[pbc3$bstat != "censor"]))
# At event time 1 (182), all 9 are at risk, (3,3,2,1) in initial states 1-4
atrisk <- pbc3$tstart < etime[1] & pbc3$tstop >= etime[1]  # all 9 at risk
table(pbc3$bili4[atrisk])

#  One event occurs at 182, a 2:1 transition  (1-2 to normal)
#  Risk scores for the non-death transitions are all exp(0) =1,
#  so the hazard matrix H will have second row of (1/3, -1/3, 0,0,0) and all
#  other rows are 0. 
with(subset(pbc3, tstop== 182), table(istate= bili4, state=bstat))

# The next four events are from 3:4, 3:2, 2:3, and 1:2, so also have 
#  simple transtions, i.e., no covariates so all risk scores are exp(0) =1
#
hmat <- array(0, dim=c(5,5,6))  # first 6 hazard matrices, start with 3,3,2,1
hmat[2,1,1] <- 1/3; hmat[2,2,1] <- -1/3   # new count= 4,2,2,1
hmat[3,4,2] <- 1/2; hmat[3,3,2] <- -1/2   # new count= 4,2,1,2
hmat[3,2,3] <- 1  ; hmat[3,3,3] <- -1     # new count= 4,3,0,2
hmat[2,3,4] <- 1/3; hmat[2,2,4] <- -1/3   # new count= 4,2,1,2
hmat[1,2,5] <- 1/4; hmat[1,1,5] <- -1/4   # new count= 3,3,1,2

# Event 6 is a transition from state 4 to death, at day 400
# For the shared hazard, the denominator is all those in states 1,2, or 4.
atrisk <- with(pbc3, tstart < etime[6] & tstop >= etime[6])
table(pbc3$bili4[atrisk]) # current states just before time 6
 
# The subject in state 2-4 is not considered to be at risk for a death.  
# The coxph routine assumes that the set of transitions that CAN happen = the 
# set that did happen at least once.
adata <- subset(pbc3, atrisk & bili4 != '2-4')
eta <- with(adata, .05*(age-50) + .6*(bili4=="1-2") + 1.1*(bili4 == ">4"))
cbind(adata[,c('id', 'age', 'tstop', 'bili4', 'bstat')], eta, risk=exp(eta)) 
basehaz <- 1/sum(exp(eta))
hmat[1,5,6] <- basehaz;            hmat[1,1,6] <- -basehaz
hmat[2,5,6] <- basehaz * exp(.6);  hmat[2,2,6] <- -basehaz*exp(.6)
hmat[4,5,6] <- basehaz * exp(1.1); hmat[4,4,6] <- -basehaz*exp(1.1)
# double check: sum of per-subject hazards at this time point = number of
#  events at this time point 
sum(basehaz * exp(eta)) ==1

tmat <- array(0., dim= dim(hmat))  # transition matrices
pstate <- matrix((4:0)/10, nrow=1)
for (i in 1:6) {
    tmat[,,i] <- as.matrix(Matrix::expm(hmat[,,i]))
    pstate <- rbind(pstate, pstate[i,]%*% tmat[,,i])
}

dtime <- c(1, which(surv3$time %in% etime))  # skip censored rows
aeq(surv3$pstate[dtime[1:7],1,], pstate)

#
# A function to do the above "by hand" calculations, over all time points
# It is verified for the particular fit we did, but written for
#  more generality.
# fit: a multi-state fit, with shared baselines
# istate: the inital state for each row of data
# p0: starting dist for compuation
# x0: curve for this set of covariates
#  
mysurv <- function(fit, istate, p0, x0, debug=0) {
    if (!inherits(fit, 'coxphms')) stop("invalid fit")
    smap <- fit$smap
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
        beta <- beta[!phrow,, drop=FALSE]
    }
          
    # Make the values for istate and state match the 1:2, etc of the fit,
    #  i.e., the order of fit$states
    # istate and state are used in tables, using factors makes sure the result
    #  is always the right size
    nstate <- length(fit$states)
    state <-  factor(fit$y[,3], 1:nstate)  # endpoint of a transition
    if (length(istate) != nrow(fit$y)) stop ("mismatched istate")
    istate <- factor(as.character(istate), fit$states)

    # set up output
    ntran <- ncol(smap)    # number of transitions
    utime <- sort(unique(fit$y[!is.na(state), 2])) # unique event times
    ntime <- length(utime)
    tmat <- matrix(0, nstate, nstate)  # transtion matrix at this time point
    pmat <- diag(nstate)               # product of transitions
    nrisk <- matrix(0., ntime, nstate) #number at risk
    wtrisk<- matrix(0., ntime, ntran)  # weighted number per transtion 
    nevent <- matrix(0L, ntime, nstate)  # number of events of each type
    pstate <- matrix(0L, ntime, nstate)  # probability in state
    hmat <- matrix(0., nstate, nstate)   # working matrix of hazards

    # eta is a matrix of (x for subject - x0) %*% coef, one row per subject,
    #  one column per transition
    eta <- (fit$x - rep(x0, each= nrow(fit$y))) %*% beta
    rwt <- exp(eta)  # the risk weight for each obs

    t1 <- fit$y[,1]
    t2 <- fit$y[,2]
    for (i in 1:ntime) {
        atrisk <- (t1 < utime[i] &  utime[i] <= t2) # risk set at this time
        event <- which(utime[i] == t2)  # potential events, at this time
        nrisk[i,]  <- c(table(istate[atrisk]))   # number at risk in each state
        nevent[i,] <- c(table(state[event]))
        # The linear predictor and hence the number at risk is different for
        #  every transition. Also, some will not be at risk for the transition.
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
        tmat <- as.matrix(Matrix::expm(hmat))      # transtion matrix
#        if (i >= debug) browser()
        pmat <- pmat %*% tmat
        pstate[i,] <- drop(p0 %*% pmat)
    }
    list(time=utime, nrisk=nrisk, nevent=nevent, pstate=pstate,  
         wtrisk= wtrisk, P=pmat)
}

test3 <- mysurv(fit3, pbc3$bili4, p0= 4:0/10,  x0 =50)
aeq(test3$pstate, surv3$pstate[match(test3$time, surv3$time),1,])

# Now with the full data set
fit2 <- coxph(list(Surv(tstart, tstop, bstat) ~ 1,
                   c(1:4):5 ~ age / common + shared), id= id, istate=bili4,
              data=pbc2, ties='breslow', x=TRUE)
surv2 <- survfit(fit2, newdata=list(age=50), p0=c(.4, .3, .2, .1, 0))
test2 <- mysurv(fit2, pbc2$bili4, p0= 4:0/10, fit2, x0 =50)
aeq(test2$pstate, surv2$pstate[match(test2$time, surv2$time),1,])


if (FALSE){
    # for testing, make a plot
    xfun <- function(i) {
        j <- match(test2$time[i], surv2$time)
        all.equal(test2$pstate[i,], surv2$pstate[j,1,])
    }
    plot(surv2, col=1:5, lwd=2)
    matpoints(test2$time, test2$pstate, col=1:5, pch='o')
}

