library(survival)

#
# Test the current prevalence version of the CI curve
#
tdata <- data.frame(id=c(1,1,1,1, 2,2,2, 3,3, 4,4,4,4, 5, 6, 6),
                    time=c(10,20,30,40, 5, 15,25, 20, 22, 6,18,34,50,10,15,20),
                    status=c(1,1,1,1, 1,1,1, 1,0, 1,1,1,0,0,1,0),
                    event= letters[c(1,2,3,4, 2,4,3, 2,2, 3,1,2,2,1, 1,1)],
                    wt = c(2,2,2,2, 1,1,1, 3,3, 1,1,1,1, 2, 1,1))

fit <- survfit(Surv(time, status) ~1, etype=event, id=id, weight=wt, tdata)
all.equal(as.vector(1-fit$surv),
          as.vector(matrix(c(0,0,2,3.5, 4.5, 2.5, 2.5, 2.5, 2.5, 0,0,0,
                             1,1,1,0,0, 6.5, 6.5, 6.5, 0, 2.5, 2.5, 2.5,
                             0,1,1,1,0,0,0,1, 7.5, 7.5, 1,1,
                             0,0,0,1,1,1,1,0,0,0, 6.5, 6.5)/10, ncol=4)))

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


# 1-fit$surv for time i and state j = total weight at that time/state in the
#  above table, divided by 10.

