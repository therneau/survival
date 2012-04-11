library(survival)
#
# Tests with the pspline function, to verify the prediction aspects
#
options(na.action=na.exclude)
aeq <- function(x,y, ...) all.equal(as.vector(x), as.vector(y), ...)

spfit <- coxph(Surv(time, status) ~ pspline(age) + ph.ecog, lung)

spfit2 <- coxph(Surv(time, status) ~ pspline(age) + ph.ecog, lung, x=TRUE)
x2 <- model.matrix(spfit)
all.equal(spfit2$x, x2)

keep <- (lung$age < 60)
x3 <- model.matrix(spfit, data=lung[keep,])  
attr(x3, 'assign') <- NULL #subscripting loses the assign attr below
all.equal(napredict(spfit$na.action,x2)[keep,], x3)

p2 <- predict(spfit, newdata=lung[keep,])
aeq(p2, predict(spfit)[keep])


p3 <- survfit(spfit)
p4 <- survfit(spfit, newdata=lung[1:2,])
temp <- scale(x2[1:2,], center=spfit$means, scale=FALSE)%*% coef(spfit)
aeq(p3$time, p4$time)
aeq(outer(-log(p3$surv), exp(temp), '*'), -log(p4$surv))

# Check out model.frame
spfit3 <- coxph(Surv(time, status) ~ pspline(age) + sex, lung,
                model=TRUE)  #avoid the missing value
m2 <- model.frame(spfit3, data=lung[keep,])
all.equal(m2, spfit3$model[keep,])

