#
# A tiny multi-state example
#
library(survival)
mtest <- data.frame(id= c(1, 1, 1,  2,  3,  4, 4, 4,  5, 5),
                    t1= c(0, 4, 9,  0,  2,  0, 2, 8,  1, 3),
                    t2= c(4, 9, 10, 5,  9,  2, 8, 9,  3, 11),
                    st= c(1, 2,  1, 2,  3,  1, 3, 0,  2,  0))

mtest$state <- factor(mtest$st, 0:3, c("censor", "a", "b", "c"))
mtest <- mtest[c(1,3,2,4,5,7,6,10, 9, 8),]  #not in time order

mfit <- survfit(Surv(t1, t2, state) ~ 1, mtest, id=id)

# True results
#
#time       state                    probabilities
#         entry  a   b  c         entry  a    b     c
#
#0        124                      1     0    0     0
#1+       1245
#2+       1235   4                3/4   1/4   0     0    4 -> a, add 3
#3+       123    4   5            9/16  1/4  3/16   0    5 -> b
#4+        23    14  5            6/16  7/16 3/16   0    1 -> a
#5+        3     14  5            3/16  7/16 6/16   0    2 -> b, exits
#8+        3     1   5  4         3/16  7/32 6/16  7/32  4 -> c
#9+                  15            0     0  19/32 13/32  1->b, 3->c & exit
# 10+            1   5                19/64 19/64 13/32  1->a

# In mfit, the "entry" state is last in the matrices
all.equal(mfit$n.risk, matrix(c(0,1,1,2,2,1,0,0,
                                0,0,1,1,1,1,2,1,
                                0,0,0,0,0,1,0,0,
                                4,4,3,2,1,1,0,0), ncol=4))
all.equal(mfit$prev,  matrix(c(8,  8, 14, 14, 7, 0,  9.5, 9.5, 
                                0,  6,  6, 12, 12,19,9.5, 9.5, 
                                0,  0,  0,  0, 7, 13, 13, 13,
                               24, 18, 12,  6, 6, 0, 0,  0)/32, ncol=4))
all.equal(mfit$n.event, matrix(c(1,0,1,0,0,0,1,0,
                                 0,1,0,1,0,1,0,0,
                                 0,0,0,0,1,1,0,0,
                                 0,0,0,0,0,0,0,0), ncol=4))
all.equal(mfit$time, c(2, 3, 4, 5, 8, 9, 10, 11))

