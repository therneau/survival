#
# This simple one is created by the __ package.  It's a sneaky way to get the
#  same variable on both sides.
# The loglik does not converge before exp(eta) overflows.
#

ysurv <- Surv(colon$time, colon$status)
ctest <- coxph(ysurv ~ status, colon)
