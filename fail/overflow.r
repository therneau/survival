# This data set has 240 thousand lines and 145 thousand events, but only
#  48 unique event times.
# Variable v2 is terribly skewed, the first iteration oversteps, and the 
#  weights on the second iteration lead to failure of the imat
#  calculation.


adata <- read.csv('overflow.csv')

test0 <- coxph(Surv(start, end, status) ~ v2, adata, iter=0)
dt <- coxph.detail( coxph(Surv(start, end, status) ~ v2,
                          adata, iter=0))
beta1 <- sum(dt$score)/sum(dt$imat)

test1 <- coxph(Surv(start, end, status) ~ v2, adata, iter=1)

test3 <- coxph(Surv(start, end, status) ~ v3, adata)

# 0 iterations works
test4 <- coxph(Surv(start, end, status) ~ v2, adata, iter=0)


# Iteration by hand
dt <- coxph.detail( coxph(Surv(start, end, status) ~ v2x,
                          adata, iter=0))
beta1 <- sum(dt$score)/sum(dt$imat)
test6 <- coxph(Surv(start, end, status) ~ v2x, init=beta1, iter=0,
               adata)

