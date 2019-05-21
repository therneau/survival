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
 
    if (missing(risk)) risk <- rep(1.0, n)
    else if (length(risk) != n || any(risk <=0)) stop("invalid risk")

    if (missing(wt)) wt <- rep(1.0, n)
    else if (length(wt) !=n || any(wt <=0)) stop("invalid wt")

    event <- y[,3] > 0
    dtime <- sort(unique(y[event,2]))  # all the event times
  
    # get the current state, and the list of transtions
    #  transitions to censor don't count
    mcheck <- multicheck(y~1, id= id)
    istate <- mcheck$istate
    temp <- attr(y, 'states')[y[event,3]]
    tmat <- table(y[event,2], from=istate[event], to=temp)
    tmat2 <- tapply(wt[event], list(y[event,2], from=mcheck$istate[event],
                                    to=temp), sum)
     
    # the number at risk in each state
    ntime <- length(dtime)
    nrisk <- wtrisk <- matrix(0, ntime, dim(tmat2)[2])
    for (i in 1:ntime) {
        atrisk <- ((y[,1] < dtime[i]) & (y[,2] >= dtime[i]))
        nrisk[i,] <- table(istate[atrisk])
        tmp <- tapply((risk*wt)[atrisk], istate[atrisk], sum)
        wtrisk[i,] <- ifelse(is.na(tmp), 0, tmp)
    }

    haz <- tmat2/c(wtrisk)
    haz <- ifelse(is.na(haz), 0, haz)
    chaz<- apply(haz, 2:3, cumsum)
    
    # make it two dimensional
    keep <- which(chaz[ntime,,] > 0)  #non-zero hazards
    klab <- outer(dimnames(haz)[[2]], dimnames(haz)[[3]], paste, sep=':')
    haz2 <- matrix(haz, nrow=ntime, dimnames=list(NULL, klab))[,keep]
    chaz2<- matrix(chaz, nrow=ntime, dimnames=list(NULL, klab))[,keep]

    list(time=dtime, nrisk=nrisk, nevent=apply(tmat, c(1,3),sum),
         haz=haz2, cumhaz=chaz2, states=mcheck$states)
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

# a coefficient of log(2)
mfit <- coxph(Surv(t1, t2, state) ~x, iter=0, mtest, id=id,
              init=rep(log(2), 6))
    
msurv <- survfit(mfit, newdata=list(x=0))

if (FALSE) {
    # this graph is very useful
    temp <- multicheck(Surv(t1, t2, state) ~1, mtest, id=id)
    plot(c(0,11), c(1,5.1), type='n', xlab="Time", ylab= "Subject")
    with(mtest, segments(t1+.1, id, t2, id, col=as.numeric(temp$istate)))
    event <- subset(mtest, state!='censor')
    text(event$t2, event$id+.2, as.character(event$state))
}
         
