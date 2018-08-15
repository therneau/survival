#
# Example from Michael Tsagris, using the data from
# Rosenwald A, Wright G, Wiestner A, Chan WC, Connors JM, Campo E, et al. 
#  The proliferation gene expression signature is a quantitative integrator
#  of oncogenic events that predicts survival in mantle cell lymphoma. 
#  Cancer cell 2003;3:185– 197.
#
# There are 92 observations, 28 events, and 8810 predictors.  Numbers 8328 and
#   8674 killed survreg.
#
load('rosenwald.rda')
library(survival)

fit0 <- survreg(Surv(days, status) ~ 1)
fit1 <- survreg(Surv(days, status) ~ x[, 8328:8329])  # worked okay
fit2 <- survreg(Surv(days, status) ~ x[, 8328])  #  not okay
fit3 <- survreg(Surv(days, status) ~ x[, 8328], init=c(7, -.21))  # okay
fit4 <- survreg(Surv(days, status) ~ x[, 8328], init=c(7, 0))  # okay

# This is a case of a horrible step in survreg's iteration path, due to a
#  bad intial guess.  

# He also sent a coxph failure

toss <- c(1, 15,22, 30, 36, 57, 63, 75, 89)
ivar <- c(4096, 6474, 1268, 39, 458, 7620, 7606, 4434, 4667, 2916, 7103, 324, 
          2795, 1594)
cfit1 <- coxph(Surv(days, status) ~ x[,ivar], subset=(-toss))
cfit2 <- coxph(Surv(days, status) ~ x[,c(7,ivar)], subset=(-toss)) #not okay

# Digging in, at iteration 9 the liner predictor has a range of -342 to 465
# The coefficients at iteration 10 leads to multiple NA elements in
#  the information matrix, but do not quite force the loglik to NA.
# The code now gives up if either the loglik OR any elements of
#  imat are non-finite.
