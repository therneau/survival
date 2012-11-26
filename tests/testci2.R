library(survival)

#
# Test the multi-state version of the CI curve
#
tdata <- data.frame(id=c(1,1,1,1, 2,2,2, 3,3, 4,4,4,4, 5, 6, 6),
                    time1=c(0, 10,20,30, 0, 5, 15,  0, 20, 0, 6,18,34, 0, 0,15),
                    time2=c(10,20,30,40, 5, 15,25, 20, 22, 6,18,34,50,10,15,20),
                    status=c(1,1,1,1, 1,1,1, 1,0, 1,1,1,0,0,1,0),
                    event= letters[c(1,2,3,4, 2,4,3, 2,2, 3,1,2,2,1, 1,1)],
                    wt = c(2,2,2,2, 1,1,1, 3,3, 1,1,1,1, 2, 1,1))
tdata$stat2 <- factor(tdata$status * as.numeric(tdata$event),
                      labels=c(" ", levels(tdata$event)))

if (FALSE) {
    # a plot of the data that helped during creation of the example
    plot(c(0,50), c(1,6), type='n', xlab='time', ylab='subject')
    with(tdata, segments(time1, id, time2, id))
    with(tdata, text(time2, id, as.numeric(stat2)-1, cex=1.5, col=2))
}
         
fit <- survfit(Surv(time1, time2, stat2) ~1, id=id, weight=wt, tdata)
true <- matrix(c(0,0,2,3.5, 4.5, 2.5, 2.5, 2.5, 2.5, 0,0,0,
                 1,1,1,0,0, 6.5, 6.5, 6.5, 0, 2.5, 2.5, 2.5,
                 0,1,1,1,0,0,0,1, 7.5, 7.5, 1,1,
                 0,0,0,1,1,1,1,0,0,0, 6.5, 6.5)/10, ncol=4)

all.equal(as.vector(fit$prev), as.vector(true))


# The exact figures for testci2.

# The subject data of  id, weight, (transition time, transition)

#1: 2 (10, 0->1)  (20, 1->2)  (30, 2->3)  (40, 3->4)
#2: 1 ( 5, 0->2)  (15, 2->4)  (25, 4->3)
#3: 3 (20, 0->2)  (22, censor)
#4: 1 ( 6, 0->3)  (18, 3->1)  (34, 1->2)  (50, censor)
#5: 2 (10, censor)
#6: 1 (15, 0->1)  (20, censor)

#Each line below follows a subject through time as a (state, weight) pair

#        |   Time
#Subject |  5     6   10   10+     15    18      20      20+    22+    25    
#----------------------------------------------------------------------
#1         0,2   0,2  1,2  1,2    1,2    1,2    2,2      2,2    2,6.5  2,6.5
#2         2,1   2,1  2,1  2,1    4,1    4,1    4,1      4,1    4,1    3,1
#3         0,3   0,3  0,3  0,4.5  0,4.5  0,4.5  2,4.5    2,4.5  
#4         0,1   3,1  3,1  3,1    3,1    1,1    1,1      1,2.5  1,2.5  1,2.5
#5         0,2   0,2  0,2  
#6         0,1   0,1  0,1  0,1.5  1,1.5  1,1.5  1,1.5    
#
# 
#          30      34    40      50
#--------------------------------------
#1         3,6.5  3,6.5  4,6.5   4,6.5
#2         3,1    3,1    3,1     3,1
#4         1,2,5  2,2.5  2,2.5 


# fit$prev for time i and state j = total weight at that time/state in the
#  above table, divided by 10.

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

