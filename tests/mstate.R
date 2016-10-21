#
# A tiny multi-state example
#
library(survival)
aeq <- function(x,y) all.equal(as.vector(x), as.vector(y))
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


# Somewhat more complex.
#  Scramble the input data
#  Not everyone starts at the same time or in the same state
#  Two "istates" that vary, only the first should be noticed.
#
tdata <- data.frame(id= c(1, 1, 1,  2,  3,  4, 4, 4,  5,  5),
                    t1= c(0, 4, 9,  1,  2,  0, 2, 8,  1,  3),
                    t2= c(4, 9, 10, 5,  9,  2, 8, 9,  3, 11),
                    st= c(1, 2,  1, 2,  3,  1, 3, 0,  3,  0),
                    i0= c(4, 4,  4, 1,  4,  4, 4, 1,  2,  2))

tdata$st <- factor(tdata$st, c(0:4),
                    labels=c("censor", "1", "2", "3", "entry"))

tfun <- function(wt, data=tdata) {
    reorder <- c(10, 9, 1, 2, 5, 4, 3, 7, 8, 6)
    new <- data[reorder,]
    new$wt <- rep(wt,length=10)[reorder]
    new
}

# These weight vectors are in the order of tdata
# w[9] is the weight for subject 5 at time 1.5, for instance
p0 <- function(w) c(w[4], w[9], 0, w[1]+ w[6])/ (w[1]+ w[4] + w[6] + w[9])

#  aj2 = Aalen-Johansen matrix at time 2, etc.
aj2 <- function(w) {
    rbind(c(1, 0, 0, 0),    # state a (1) stays put
          c(0, 1, 0, 0),
          c(0, 0, 1, 0),
          c(w[6], 0, 0, w[1])/(w[1] + w[6]))  #subject 4 moves to 'a'
}
aj3 <- function(w) rbind(c(1, 0, 0, 0),   
                         c(0, 0, 1, 0),  # 5 moves from b to c
                         c(0, 0, 1, 0),
                         c(0, 0, 0, 1))
aj4 <- function(w) rbind(c(1, 0, 0, 0),
                         c(0, 1, 0, 0),  
                         c(0, 0, 1, 0),
                         c(w[1], 0, 0, w[5])/(w[1] + w[5])) #1 moves from 4 to a
aj5 <- function(w) rbind(c(w[2]+w[7], w[4], 0, 0)/(w[2]+ w[4] + w[7]), #2 to b
                         c(0, 1, 0, 0),  
                         c(0, 0, 1, 0),
                         c(0, 0, 0, 1))
aj8 <- function(w) rbind(c(w[2], 0, w[7], 0)/(w[2]+ w[7]), # 4  to c
                         c(0, 1, 0, 0),  
                         c(0, 0, 1, 0),
                         c(0, 0, 0, 1))
aj9 <- function(w) rbind(c(0, 1, 0, 0), # 1  to b
                         c(0, 1, 0, 0),  
                         c(0, 0, 1, 0),
                         c(0, 0, 1 ,0)) # 3 to c
aj10 <- function(w)rbind(c(1, 0, 0, 0),
                         c(1, 0, 0, 0),  #1 back to a
                         c(0, 0, 1, 0),
                         c(0, 0, 0, 1))

#time       state               
#         a   b  c  entry
#
#1        2   5     14       initial distribution
#2        24  5     1        4 -> a, add 3
#3        24     5  13       5 from b to c
#4       124     5   3       1 -> a
#5        14     5   3       2 -> b, exits
#8        1      45  3       4 -> c
#9            1  45          1->b, 3->c & exit
#10       1      45          1->a

# The prevalence is a product of matrices
doprev <- function(w) {
    p1 <- p0(w)
    p2 <- p1 %*% aj2(w)
    p3 <- p2 %*% aj3(w)
    p4 <- p3 %*% aj4(w)
    p5 <- p4 %*% aj5(w)
    p8 <- p5 %*% aj8(w)
    p9 <- p8 %*% aj9(w)
    p10<- p9 %*% aj10(w)
    rbind(p2, p3, p4, p5, p8, p9, p10, p10)
}

# Check the pstate estimate
w1 <- rep(1, 10)
mtest2 <- tfun(w1)
mfit2 <- survfit(Surv(t1, t2, st) ~ 1, tdata, id=id, istate=i0) #rdered
all.equal(mfit2$prev, doprev(w1))
aeq(mfit2$p0, p0(w1))

mfit2b <- survfit(Surv(t1, t2, st) ~ 1, mtest2, id=id, istate=i0)#scrambled
all.equal(mfit2b$prev, doprev(w1))
aeq(mfit2b$p0, p0(w1))

mfit2b$call <- mfit2$call <- NULL
all.equal(mfit2b, mfit2) 

# The derivative of a matrix product AB is (dA)B + A(dB) where dA is the
#  elementwise derivative of A and etc for B.
# dp0 creates the derivatives of p0 with respect to each subject, a 5 by 4
#  matrix
dp0 <- function(w) {
  p <- p0(w)
  w0 <- w[c(1,4,6,9)]  # the 4 obs at the start, subjects 1, 2, 4, 5
  rbind(c(0, 0, 0, 1) - p,   # subject 1 affects p[4]
        c(1, 0, 0, 0) - p,   # subject 2 affects p0[1]
        0,                   # subject 3 affects none
        c(0, 0, 0, 1) - p,   # subject 4 affect p[4]
        c(0, 1, 0, 0) - p) / sum(w0)
}
  
mfit3 <- survfit(Surv(t1, t2, st) ~ 1, tdata, id=id, istate=i0,
                 influence=TRUE)
