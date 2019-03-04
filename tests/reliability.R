options(na.action=na.exclude) # preserve missings
options(contrasts=c('contr.treatment', 'contr.poly')) #ensure constrast type
aeq <- function(x,y, ...) all.equal(as.vector(x), as.vector(y), ...)

library(survival)
# Many of these results are found in Meeker and Escobar, Statistical
#  Methods for Reliability Data, 1998.

# Generator fan- the simplest possible data
# Meeker & Escobar, page 166 claims theta=10.26476, but our solution of
#  10.26477 is a touch better.
g1 <- survreg(Surv(hours, status) ~1, genfan, dist="exponential")
aeq(round(g1$coef,5), 10.26477)  
g2 <- survreg(Surv(hours, status) ~1, genfan, dist="exponential",
              init=10.26476, iter=0)
g1$loglik - g2$loglik  # larger loglik by a whisker

# loglik on a transformed scale changes by the Jacobian of the 
#  transform, but only for non-censored observations
g3 <- survreg(Surv(log(hours), status) ~1, genfan,
              dist="extreme", scale=1)
aeq(g1$coef, g3$coef)
aeq(g1$loglik, g3$loglik - sum(log(genfan$hours[genfan$status==1])))

# verify the loglik
lambda <- exp(-g1$coef)
elog <- with(genfan, ifelse(status==1, dexp(hours, lambda, TRUE), 
                                       pexp(hours, lambda, FALSE, TRUE)))
aeq(g1$loglik[1], sum(elog))


# Glass capacitor data: 2 covariates
nr <- nrow(capacitor)
table17 <- with(capacitor, tapply(1:nr, list(temperature, voltage), 
              function(x) {
                  tfit <- survreg(Surv(time, status)~1, dist="weibull", 
                                  data=capacitor, subset=x)
                  tfit$coef[1]
              }))

# the mu-hat elements of Table 17.3
meeker17 <- matrix(c(7.13, 7.01, 7.10, 6.28, 6.57, 6.00, 6.54, 6.25), nrow=2)
aeq(round(table17, 2), meeker17)

# This will match table 17.4
fitc1 <- survreg(Surv(time, status)~ voltage + temperature, 
                 dist="weibull", capacitor)
fitc2 <- update(fitc1, . ~ . + voltage * temperature)

cap2 <- capacitor
cap2$voltage <- cap2$voltage-250
cap2$temperature <- cap2$temperature - 170
x1 <- survreg(Surv(time, status)~ voltage + temperature, 
                 dist="weibull", cap2)
x2 <- update(x1, . ~ . + voltage * temperature) 

# other models
fitig <- survreg(Surv(fail)~ voltage + temperature,  
	dist = "gaussian", data = capacitor)
summary(fitig)

fitix <- survreg(Surv(fail)~ voltage + temperature,  
	dist = "extreme", data = capacitor)
summary(fitix)

fitil <- survreg(Surv(fail)~ voltage + temperature,  
	dist = "logistic", data = capacitor)
summary(fitil)
