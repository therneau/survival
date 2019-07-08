library(survival)
aeq <- function(x, y) all.equal(as.vector(x), as.vector(y))
# 
# Compute the hazard functions for a multi-state Cox model
#
coxhaz <- function(y, id, risk, wt) {
    # y should be a multi-state survival
    if (!inherits(y, "Surv") || attr(y, "type") != "mcounting")
        stop("invalid response")

    n <- nrow(y)
    if (missing(id) || length(id) !=n) stop("invalid id")
     
    if (missing(wt)) wt <- rep(1.0, n)
    else if (length(wt) !=n || any(wt <=0)) stop("invalid wt")

    # get the current state, and the list of transtions
    #  transitions to censor don't count
    mcheck <- survcheck(y~1, id= id)
    states <- mcheck$states
    nstate <- length(states)
    istate <- mcheck$istate

    event <- y[,3] > 0
    temp <- attr(y, 'states')[y[event,3]]
    tmat <- table(y[event,2], from=istate[event], to=temp)
    tmat2 <- tapply(wt[event], list(y[event,2], from=mcheck$istate[event],
                                    to=temp), sum)
    tmat2 <- ifelse(is.na(tmat2), 0, tmat2)


    # Hazards can be done one at a time.  For each of them the risk
    #  weight vector for the subjects can be different.
    # First organized the material as a 2 dim matrix
    temp <- apply(tmat, 2:3, sum)
    keep <- which(temp>0)
    from <- states[row(temp)[keep]]
    hlab <- outer(rownames(temp), colnames(temp), paste, sep=':')[keep]
    nhaz <- length(keep)
    nevent <- matrix(tmat2, nrow(tmat2))[,keep]
    dtime <- sort(unique(y[event,2]))
    ntime <- length(dtime)

    if (missing(risk)) risk <- matrix(1, nrow=n, ncol=nhaz)
    if (!is.matrix(risk) || nrow(risk) != n || ncol(risk) != nhaz)
        stop("invalid risk matrix")
    risk <- risk * wt

    # get the weighted at risk at each time
    wtrisk <- matrix(, length(dtime), nhaz)
    statematch <- outer(istate, from, function(x, y) x==y)
    risk <- ifelse(statematch, risk, 0)
    for (i in 1:ntime) {
        atrisk <- (y[,1]< dtime[i] & y[,2] >= dtime[i])
        wtrisk[i,] <- colSums(risk[atrisk,, drop=FALSE]) 
    }

    haz <- nevent/ifelse(wtrisk==0, 1, wtrisk)   # avoid 0/0
    chaz<- apply(haz, 2, cumsum)
    
    list(time=dtime, nrisk=wtrisk, nevent=nevent,
         haz=haz, cumhaz=chaz, states=states)
}

mtest <- data.frame(id= c(1, 1, 1,  2,  3,  4, 4, 4,  5, 5),
                    t1= c(0, 4, 9,  0,  2,  0, 2, 8,  1, 3),
                    t2= c(4, 9, 10, 5,  9,  2, 8, 9,  3, 11),
                    state= c(1, 2,  1, 2,  3,  1, 3, 0,  2,  0),
                    x = c(0, 0,  0, 1,  1,  0, 0, 0,  2,  2))
mtest$state <- factor(mtest$state, 0:3, c("censor", "a", "b", "c"))

# True results
#
#time       at risk               events
#         entry  a   b  c        
#
#2        1245                   4 -> a
#3        1235   4               5 -> b
#4        123    4   5           1 -> a
#5         23    14  5           2 -> b, exits
#8         3     14  5           4 -> c
#9         3     1   5  4        1->b, 3->c & exit, 4 censored
#10                  15          1->a, exit
#11                   5          censor

# with all coefficients =0 
check1 <- coxhaz(Surv(mtest$t1, mtest$t2, mtest$state), mtest$id)
fit1 <- survfit(Surv(t1, t2, state) ~1, mtest, id=id)
aeq(check1$cumhaz, fit1$cumhaz[match(check1$time, fit1$time),])

dummy <- data.frame(x=1:2)
cox0 <-  coxph(Surv(t1, t2, state) ~x, iter=0, mtest, id=id)
cfit0 <- survfit(cox0, newdata=dummy)
indx <- match(check1$time, cfit0$time)
aeq(check1$cumhaz, cfit0$cumhaz[indx,1,])
aeq(check1$cumhaz, cfit0$cumhaz[indx,2,])

# a fixed coefficient
mfit <- coxph(Surv(t1, t2, state) ~x, iter=0, mtest, id=id,
              init= log(1:6))
msurv <- survfit(mfit, newdata=list(x=0:1))
mrisk <- exp(outer(mtest$x, log(1:6), '*'))  # hazards for each transition
check2 <- with(mtest, coxhaz(Surv(t1, t2, state), id=id, risk=mrisk))
aeq(check2$cumhaz, msurv$cumhaz[indx,1,])

# a different predicted x multiplies each column of the cumulative hazard
aeq(check2$cumhaz %*% diag(1:6), msurv$cumhaz[indx,2,]) 

if (FALSE) {
    # this graph is very useful
    temp <- survcheck(Surv(t1, t2, state) ~1, mtest, id=id)
    plot(c(0,11), c(1,5.1), type='n', xlab="Time", ylab= "Subject")
    with(mtest, segments(t1+.1, id, t2, id, col=as.numeric(temp$istate)))
    event <- subset(mtest, state!='censor')
    text(event$t2, event$id+.2, as.character(event$state))
}
         

# slight change, add a few censored subjects
#  all the events happen on even numbered days
test2 <- data.frame(id= c(1, 1, 1,  2,  3,  4, 4, 4,  5, 5,
                          6, 7, 8,  9),
                    t1= c(0, 8, 18,  0,  4,  0, 4, 16,  2, 6,
                          0, 0, 7,  8),
                    t2= c(8, 18, 20, 10,  18,  4, 16, 18,  6, 22,
                          5, 10,  10, 15),
                    state= c(1, 2,  1, 2,  3,  1, 3, 0,  2,  0,0,0,0,0),
                    x = c(0, 0,  0, 1,  1,  0, 0, 0,  2,  2, 1, 1, 2, 0))
test2$state <- factor(test2$state, 0:3, c("censor", "a", "b", "c"))

cox2 <- coxph(Surv(t1,t2, state) ~ x, id=id, test2, 
                     init=log(1:6), iter=0)
csurv2 <- survfit(cox2, newdata=data.frame(x=0:1))
