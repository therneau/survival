library(survival)

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
                      labels=c(" ", levels(tdata$event)))

if (FALSE) {
    # a plot of the data that helped during creation of the example
    plot(c(0,50), c(1,6), type='n', xlab='time', ylab='subject')
    with(tdata, segments(time1, id, time2, id))
    with(tdata, text(time2, id, as.numeric(stat2)-1, cex=1.5, col=2))
}
         
fit <- survfit(Surv(time1, time2, stat2) ~1, id=id, weight=wt, tdata)

truth <- matrix(c(0,0,2, 3.5, 4.5, rep(2.25, 5),   0, 0,
                  1,1,1, 0,   0,   6.75 ,6.75, 6.75,  0, 9/4, 9/4 , 9/4,
                  0,1,1, 1,   1,   0, 0,    1,    7.75, 7.75, 0, 0, 
                  0,0,0, 0,   1,   1, 1,    0,     0,    0,   7.75, 7.75)/10,
                ncol=4)
#dimnames(truth) <- list(c(5,6, 10, 15, 18, 20, 25, 30, 34, 40, 50), 
#                        c('a', 'b', 'c', 'd'))
all.equal(fit$prev, truth)


# The exact figures for testci2.

# The subject data of  id, weight, (transition time, transition)

#1: 2 (10, 0->1)  (20, 1->2)  (30, 2->3)  (40, 3->4)  no data after 40=censored
#2: 1 ( 5, 0->2)  (15, 2->4)  (25, 4->3)  no data after 25 implies censored then
#3: 3 (20, 0->2)  (22, censor)
#4: 1 ( 6, 0->3)  (18, 3->1)  (34, 1->2)  (50, censor)
#5: 2 (10, censor)
#6: 1 (15, 0->1)  (20, censor)

# Each line below follows a subject through time as a (state, rdist weight) pair
#  using the redistribute to the right algorithm.
# RDR algorithm: at each censoring (or last fu) a subject's weight is put into
#  a "pool" for that state and their weight goes to zero. The pool is
#  dynamically shared between all members of the state proportional to their
#  original case weight.  If someone new enters the state it is reapportioned,
#  when someone leaves they take their current portion to the new state with them.
# When there is no one left in the state the last one to be censored carries
#  it forward as a "ghost", perhaps temporarily, these are the () entries below.
#
#Subject |   Time
#  t     |  5     6   10   10+     15      18      20       20+       22+   
#---------------------------------------------------------------------------
#1  2      -,0   -,0  a,0  a,0    a,1/3    a,1/4  b,14/20   b,7/10    b,19/4
#2  1      b,0   b,0  b,0  b,0    d,0      d,0    d,0       d,0       d,0   
#3  3      -,0   -,0  -,0  -,3/2  -,3/2    -,3/2  b,21/20   b,21/20
#4  1      -,0   c,0  c,0  c,0    c,0      a,1/8  a,1/8     a,5/4     a,5/4
#5  2      -,0   -,0  -,0  
#6  1      -,0   -,0  -,0  -,1/2  a,1/6    a,1/8  a,1/8    
#
# 
#          25      25+     30        34      40       50
#--------------------------------------------------------
#1        b,19/4   b,19/4  c,23/4   c,23/4   d,23/4  (d, 23/4)
#2        c,0      (c, 1)
#4        a,5/4    a,5/4    a,5/4   b,5/4    b,5/4   b,5/4
#
#
# fit$prev for time i and state j = total weight at that time/state in the
#  above table (original weight + redistrib), divided by 10.

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

