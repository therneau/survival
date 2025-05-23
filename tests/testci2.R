library(survival)
aeq <- function(x,y, ...) all.equal(as.vector(x), as.vector(y), ...)
#
# Test the multi-state version of the CI curve
#
tdata <- data.frame(id=c(1,1,1,1, 2,2,2, 3,3, 4,4,4,4, 5, 6, 6),
                    time1=c(0, 10,20,30, 0, 5, 15,  0, 20, 0, 6,18,34, 0, 0,15),
                    time2=c(10,20,30,40, 5, 15,25, 20, 22, 6,18,34,50,10,15,20),
                    status=c(1,1,1,1, 1,1,1, 1,0, 1,1,1,0,0,1,0),
                    event= letters[c(1,2,3,4, 2,4,3, 2,2, 3,1,2,2,1, 1,1)],
                    wt = c(2,2,2,2, 1,1,1, 3,3, 1,1,1,1, 2, 1,1),
                    stringsAsFactors=TRUE)
tdata$stat2 <- factor(tdata$status * as.numeric(tdata$event),
                      labels=c("censor", levels(tdata$event)))
         
fit <- survfit(Surv(time1, time2, stat2) ~1, id=id, weights=wt, tdata,
               influence=TRUE)

# The exact figures for testci2.
# The subject data of  id, weight, (transition time, transition)

#1: 2 (10, 0->a)  (20, a->b)  (30, b->c)  (40, c->d)  no data after 40=censored
#2: 1 ( 5, 0->b)  (15, b->d)  (25, d->c)  no data after 25 implies censored then
#3: 3 (20, 0->b)  (22, censor)
#4: 1 ( 6, 0->c)  (18, c->a)  (34, a->b)  (50, censor)
#5: 2 (10, censor)
#6: 1 (15, 0->a)  (20, censor)

# Each line below follows a subject through time as a (state, rdist weight) pair
#  using the redistribute to the right algorithm.
# RDR algorithm: at each censoring (or last fu) a subject's weight is put into
#  a "pool" for that state and their weight goes to zero. The pool is
#  dynamically shared between all members of the state proportional to their
#  original case weight, when someone leaves they take their portion of the
#  pool to the new state.

# Table of case weights and state, blank is weight of zero
#    time   5    6   10   15    18    20     25    30    34    40    50
# -----------------------------------------------------------------------
# id, wt
#  1, 2     -    -    a    a      a    b      b     c     c     d    
#  2, 1     b    b    b    d      d    d      c    
#  3, 3     -    -    -    -      -    b      
#  4, 1     -    c    c    c      a    a      a     a     b     b     b
#  5, 2     -    -    -
#  6, 1     -    -    -    a      a    a      

# Pool weights
#         10  10+   15    18    20    20+   22+   25   25+   30   34    40   40+
#    -     0   2    3/2  3/2    0 
#    a     0   0    1/2  1/2    1/4   5/4   5/4   5/4  5/4   5/4
#    b     0   0     0    0     7/4   7/4  19/4  19/4 19/4        5/4   5/4  5/4
#    c     0   0     0    0     0                      1    23/4 23/4
#    d     0   0     0    0     0                                      23/4 31/4

# fit$pstate for time i and state j = total weight at that time/state in the
#  above table (original weight + redistrib), divided by 10.

# time            5  6   10    15    18    20     25    30    34    40    50
truth <- matrix(c(0, 0,   2,    3,    4,    2,     1,    1,    0,    0,    0,
                  1, 1,   1,    0,    0,    5,     2,    0,    1,    1,    1,
                  0, 1,   1,    1,    0,    0,     1,    2,    2,    0,    0,
                  0, 0,   0,    1,    1,    1,     0,    0,    0,    2,    0) +
                c(0, 0,   0,   .5,   .5,   1/4,   5/4,  5/4,   0,    0,    0,
                  0, 0,   0,    0,    0,   7/4,  19/4,   0,   5/4,   5/4,  5/4,
                  0, 0,   0,    0,    0,    0,     0,  23/4, 23/4,   0,    0, 
                  0, 0,   0,    0,    0,    0,     0,    0,    0,  23/4,  31/4),
                ncol=4)
truth <- truth[c(1:6, 6:11),]/10  #the explicit censor at 22

#dimnames(truth) <- list(c(5, 6, 10, 15, 18, 20, 25, 30, 34, 40, 50),
#                        c('a', 'b', 'c', 'd')
aeq(truth, fit$pstate[,2:5])

# Test the dfbetas
# It was a big surprise, but the epsilon where a finite difference approx to
#  the derivative is most accurate is around 1e-7 = approx sqrt(precision).
# Smaller eps makes the approximate derivative worse.
# There is a now a formal test in mstate.R, not approximate.

# compute the per observation influence first
n <- nrow(tdata)     
U <- array(0, dim=c(n, dim(fit$pstate)))
eps <- sqrt(.Machine$double.eps)
n <- nrow(tdata)     
for (i in 1:n) {
    twt <- tdata$wt
    twt[i] <- twt[i] + eps
    tfit <- survfit(Surv(time1, time2, stat2) ~ 1, id=id, tdata,
                    weights=twt)
    U[i,,] <- (tfit$pstate - fit$pstate)/eps  #finite difference approx
}
dfbeta <- rowsum(tdata$wt*matrix(U,nrow=n), tdata$id) # per subject
dfbeta <- array(dfbeta, dim=c(6,12,5))
aeq(dfbeta, fit$influence, tolerance= eps*10)

aeq(fit$std.err, sqrt(apply(fit$influence.pstate^2, 2:3, sum)))

if (FALSE) {
    # a plot of the data that helped during creation of the example
    plot(c(0,50), c(1,6), type='n', xlab='time', ylab='subject')
    with(tdata, segments(time1, id, time2, id))
    with(tdata, text(time2, id, as.numeric(stat2)-1, cex=1.5, col=2))
}

if (FALSE) {
    # The following lines test out 4 error messages in the routine
    #
    # Gap in follow-up time, id 2
    survfit(Surv(c(0,5,9,0,5,0), c(5,9,12, 4, 6, 3), factor(c(0,0,1,1,0,2))) ~1,
            id=c(1,1,1,2,2,3))
    # mismatched weights
    survfit(Surv(c(0,5,9,0,5,0), c(5,9,12, 5, 6, 3), factor(c(0,0,1,1,0,2))) ~1,
            id=c(1,1,1,2,2,3), weights=c(1,1,2,1,1,4))
    # in two groups at once
     survfit(Surv(c(0,5,9,0,5,0), c(5,9,12, 5, 6, 3), factor(c(0,0,1,1,0,2))) ~
               c(1,1,2,1,1,2), id=c(1,1,1,2,2,3)) 
    # state change that isn't a state change (went from 1 to 1)
        survfit(Surv(c(0,5,9,0,5,0), c(5,9,12, 5, 6, 3), factor(c(0,1,1,1,0,2))) ~1,
            id=c(1,1,1,2,2,3))
}

# Check the start.time option
#
# Later work showed this test has to be false.  At time 0 everyone starts in
# state (s0), but by time 20 many have shifted to another.  fit2 picks up at
# the right place, but because there is no istate varaible, fit2x starts
# everyone in (s0) at time 20.  There is no way for survfit to know.
if (FALSE) {
fit2 <- survfit(Surv(time1, time2, stat2) ~1, id=id, weights=wt, tdata,
                start.time=20)
data2 <- subset(tdata, time2>= 20)
fit2x <- survfit(Surv(time1, time2, stat2) ~1, id=id, weights=wt, data2)

ii <- names(fit2)[!(names(fit2) %in%  c("call", "start.time"))]
all.equal(unclass(fit2)[ii], unclass(fit2x)[ii])
}
