#
library(survival)

# Some simple tests of the Surv function
#  The first two are motivated by a bug, pointed out by Kevin Buhr,
#    where a mixture of NAs and invalid values didn't work right
#  Even for the simplest things a test case is good.
#  All but the third should produce warning messages
aeq <- function(x,y) all.equal(as.vector(x), as.vector(y))
temp <- Surv(c(1, 10, 20, 30), c(2, NA, 0, 40), c(1,1,1,1))
aeq(temp, c(1,10,NA,30,  2,NA,0,40, 1,1,1,1))

temp <- Surv(c(1, 10, 20, 30), c(2, NA, 0, 40), type='interval2')
aeq(temp, c(1,10,20,30,  2,1,1,40, 3,0,NA,3))

#No error
temp <- Surv(1:5)
aeq(temp, c(1:5, 1,1,1,1,1))

temp1 <- Surv(c(1,10,NA, 30, 30), c(1,NA,10,20, 40), type='interval2')
temp2 <- Surv(c(1,10,10,30,30), c(9, NA, 5, 20,40), c(1, 0, 2,3,3),
              type='interval')
aeq(temp1, temp2)
aeq(temp1, c(1,10,10,30,30, 1,1,1,1, 40, 1,0,2,NA,3))

# Use of inf
temp1 <- Surv(c(1,10,NA, 30, 30), c(1,NA,10,30, 40), type='interval2')
temp2 <- Surv(c(1,10,-Inf, 30, 30), c(1,Inf,10,30, 40), type='interval2')
aeq(temp1, temp2)

# Verify sorting and order routines
#  These fail in 3.4, succeed in 3.5 due to a system change in how
#  xtfrm.Surv is used.  
x1 <- Surv(c(4, 6, 3, 2, 1, NA, 2), c(1,0, NA, 0,1,1,1))
all.equal(order(x1), c(5,7, 4, 1, 2, 3, 6))
all.equal(order(x1, decreasing=TRUE), c(2,1,4,7,5, 3, 6))
all.equal(sort(x1), x1[c(5,7,4,1,2)])

x2 <- Surv(c(4, 6, 3, 2, 1, NA, 2), c(1,0, NA, 0,1,1,1), type='left')
all.equal(order(x2), c(5,4, 7, 1, 2, 3, 6))

x3 <- Surv(c(1,5,NA,7, 9), c(6, 6, 4, NA, 9), type="interval2")
all.equal(sort(x3), x3[c(1,3,2,4,5)])

x4 <- Surv(c(1,5,6,5,2, 4), c(3, 7, 7, 6, 3, NA), factor(c(1, 2, 0, 1, 1, 0)))
all.equal(sort(x4), x4[c(1, 5, 4, 2,3)])
all.equal(sort(x4, na.last=FALSE), x4[c(6,1,5,4,2,3)])
