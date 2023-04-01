# Tests of the Brier score.
#  Start with the example in the vignette
library(survival)

rott2 <- rotterdam
ignore <- with(rott2, recur ==0 & death==1 & rtime < dtime)
rott2$rfs <- with(rott2, ifelse(recur==1 | ignore, recur, death))
rott2$rfstime <- with(rott2, ifelse(recur==1 | ignore, rtime, dtime))/365.25

rsurv <- survfit(Surv(rfstime, rfs) ~1, rott2)  #KM
rfit <- coxph(Surv(rfstime, rfs) ~ pspline(age) + meno + size + pmin(nodes,12),
              rott2)

tau <- c(2,4,6, 8)  # four tau values
bfit <- brier(rfit, times=tau)

# Now by hand
wtmat <- rttright(Surv(rfstime, rfs) ~ 1, rott2, times=tau)
psurv <- survfit(rfit, newdata= rott2) # one curve per subject
yhat  <- 1- summary(psurv, times=tau)$surv
ybar  <- 1- summary(rsurv, times=tau)$surv

y <- with(rott2, cbind(rfstime <=tau[1] & rfs==1,
                       rfstime <=tau[2] & rfs==1,
                       rfstime <=tau[3] & rfs==1,
                       rfstime <=tau[4] & rfs==1)) * 1L
ss1 <- colSums(wtmat * (y - t(yhat))^2)
ss2 <- colSums(wtmat * (y - rep(ybar, each=nrow(y)))^2)

all.equal(unname(1- ss1/ss2), bfit$rsquare)
all.equal(unname(ss1), bfit$brier)
