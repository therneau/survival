#
# Make sure that useless intervals do not cause issues
#

test2 <- data.frame(time1 =c(1, 2, 5, 2, 1, 7, 3, 4, 8, 8, 3),
                    time2 =c(2, 3, 6, 7, 8, 9, 9, 9,14,17, 5),
                    event =c(1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0),
                    x     =c(1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 500) )

# The data set is the same as book3.R, except for the wild observation
#  with x=500 whose time interval of (4,5) overlaps no events.

fit1 <- coxph(Surv(time1, time2, event) ~ x, test2, subset=(x<100))
fit2 <- coxph(Surv(time1, time2, event) ~ x, test2)

ii <- match(c("coefficients", "var", "loglik", "score", "iter", 
              "wald.test", "concordance"), names(fit1))
all.equal(fit1[ii], fit2[ii])

all.equal(fit1$residuals, fit2$residuals[-11])

# The mean differs condiderably, and so to the linear predictors

# Now the same with a penalized model
fit3 <- coxph(Surv(time1, time2, event) ~ ridge(x, theta=.1), test2,
              subset= (x< 100))
fit4 <-  coxph(Surv(time1, time2, event) ~ ridge(x, theta=.1), test2)
all.equal(fit3[ii], fit4[ii])
all.equal(fit3$residuals, fit4$residuals[-11])
