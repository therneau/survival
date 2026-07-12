#
# Regression tests for two name typos in survfit methods
#
library(survival)
aeq <- function(x, y, ...) all.equal(as.vector(x), as.vector(y), ...)

#
# aggregate.survfit: the FUN= branch assigned to a stray "news" object
#  instead of "newx", raising an error for any user-supplied FUN.
#
cox1 <- coxph(Surv(time, status) ~ ph.ecog + sex, lung)
nd <- expand.grid(ph.ecog=0:2, sex=1:2)
sf <- survfit(cox1, newdata=nd)

ag1 <- aggregate(sf)                          # default FUN = mean
stopifnot(aeq(ag1$surv, rowMeans(sf$surv)))

ag2 <- aggregate(sf, FUN=median)              # explicit FUN
stopifnot(aeq(ag2$surv, apply(sf$surv, 1, median)))

#
# summary.survfit: when a newdata column name collides with an internal
#  column, the rename referenced an undefined "newdata" instead of
#  "fit$newdata", raising an error.
#
nd2 <- data.frame(ph.ecog=0:2, sex=1:2, surv=1:6)   # "surv" collides
sf2 <- survfit(cox1, newdata=nd2)
sm2 <- summary(sf2, times=c(100, 200, 300), data.frame=TRUE)
stopifnot("surv_" %in% names(sm2))            # user column renamed
stopifnot("surv"  %in% names(sm2))            # internal column retained
