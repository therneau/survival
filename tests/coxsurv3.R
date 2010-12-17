library(survival)
aeq <- function(x,y) all.equal(as.vector(x), as.vector(y))
# One more test on coxph survival curves, to test out the individual
#  option.  First fit a model with a time dependent covariate
#
test2 <- data.frame(start=c(1, 2, 5, 2, 1, 7, 3, 4, 8, 8),
                    stop =c(2, 3, 6, 7, 8, 9, 9, 9,14,17),
                    event=c(1, 1, 1, 1, 1, 1, 1, 0, 0, 0),
                    x    =c(1, 0, 0, 1, 0, 1, 1, 1, 0, 0) )

# True hazard function, from the validation document
lambda <- function(beta, x=0, method='efron') {
    r <- exp(beta)
    lambda <- c(1/(r+1), 1/(r+2), 1/(3*r +2), 1/(3*r+1),
                1/(3*r+1), 1/(3*r+2) + 1/(2*r +2))
    if (method == 'breslow') lambda[9] <- 2/(3*r +2)
    list(time=c(2,3,6,7,8,9), lambda=lambda)
    }

fit <- coxph(Surv(start, stop, event) ~x, test2)
# A curve for someone who never changes
surv1 <-survfit(fit, newdata=list(x=0), censor=FALSE)

true <- lambda(fit$coef, 0)

aeq(true$time, surv1$time)
aeq(-log(surv1$surv), cumsum(true$lambda))

# Reprise it with a time dependent subject who doesn't change
data2 <- data.frame(start=c(0, 4, 9, 11), stop=c(4, 9, 11, 17),
                      event=c(0,0,0,0), x=c(0,0,0,0))
surv2 <- survfit(fit, newdata=data2, individual=TRUE, censor=FALSE)
aeq(surv2$surv, surv1$surv)



fit <- coxph(Surv(start, stop, event) ~ age + year+ transplant, jasa1)

# A curve for someone who never changes
surv1 <- survfit(fit, newdata=list(age= -10, year=3, transplant=0))

# Reprise it with a subject who doesn't change
newdata <- data.frame(start=c(0,50,100), stop=c(50,100, max(jasa1$stop)), 
                   event=c(1,1,1), age=rep(-10,3), year=rep(3,3),
                   transplant=rep(0,3))
surv2 <- survfit(fit, newdata, individual=T)

# Have to use unclass to avoid [.survfit trying to pick curves,
#  remove the final element "call" because it won't match
all.equal(unclass(surv1)[-length(surv1)],
          unclass(surv2)[-length(surv2)])

