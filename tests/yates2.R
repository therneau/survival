library(survival)
# Tests for glm

gfit1 <- glm(skips ~ Mask* Opening + Solder, data=solder, poisson)
yg1 <- yates(gfit1, ~Mask)
yg2 <- yates(gfit1, ~ Mask, predict='response')

# Fit a model with an NA coef
# there is only one subject with ph.ecog=3, so the interaction of
#  ph.ecog=3 with sex is not estimable.
fit2 <- coxph(Surv(time, status) ~ factor(ph.ecog) *sex + age, lung)
fit2

# The marginal estimate for that group, averaged over sex, is not available
#  so becomes NA, as does the test for "all 4 ph.ecog population averages are
#  the same
yf1 <- yates(fit2, ~ph.ecog)
yf1

# Do it by hand,
temp <- subset(lung, !is.na(ph.ecog))  # the one person removed
smean <- mean(temp$sex)   # fraction males +1
amean <- mean(temp$age)
tdata <- expand.grid(ph.ecog=0:3, sex=smean, age=amean)
yhat <- predict(fit2, newdata=tdata)
all.equal(c(yhat[1:3], NA),  yf1$estimate[, "pmm"], check.attributes=FALSE)

# for predicted risk, the code uses a simulation variance
set.seed(1950)  # make the se reproducable
yf2 <- yates(fit2, ~ph.ecog, predict="risk")
yf2
# compute PMM by hand
dummy <- temp
est <- c(0,0,0, NA)
for (i in 0:2) {
    dummy$ph.ecog <- i
    est[i+1] <- mean(exp(predict(fit2, newdata=dummy)))
}
all.equal(est, yf2$estimate[,"pmm"])


# For age, we don't have missing, so NA does not appear
yates(fit2, ~ age, levels=c(50,60, 70))
