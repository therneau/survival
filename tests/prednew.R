#
# Make sure that the newdata argument works for various
#   predictions
# We purposely use a subset of the lung data that has only some
#   of the levels of the ph.ecog
library(survival)
options(na.action=na.exclude, contrasts=c('contr.treatment', 'contr.poly'))
aeq <- function(x,y) all.equal(as.vector(x), as.vector(y))

myfit <- coxph(Surv(time, status) ~ age + factor(ph.ecog) + strata(sex), lung)

keep <- which(lung$sex==1 & (lung$ph.ecog==1 | lung$ph.ecog==2))
p1 <- predict(myfit, type='lp')
p2 <- predict(myfit, type="lp", newdata=lung[keep,])
p3 <- predict(myfit, type='lp', se.fit=TRUE)
p4 <- predict(myfit, type="lp", newdata=lung[keep,], se.fit=TRUE)
aeq(p1[keep], p2)
aeq(p1, p3$fit)
aeq(p1[keep], p4$fit)
aeq(p3$se.fit[keep], p4$se.fit)

p1 <- predict(myfit, type='risk')
p2 <- predict(myfit, type="risk", newdata=lung[keep,])
p3 <- predict(myfit, type='risk', se.fit=TRUE)
p4 <- predict(myfit, type="risk", newdata=lung[keep,], se.fit=TRUE)
aeq(p1[keep], p2)
aeq(p1, p3$fit)
aeq(p1[keep], p4$fit)
aeq(p3$se.fit[keep], p4$se.fit)

p1 <- predict(myfit, type='expected')
p2 <- predict(myfit, type="expected", newdata=lung[keep,])
p3 <- predict(myfit, type='expected', se.fit=TRUE)
p4 <- predict(myfit, type="expected", newdata=lung[keep,], se.fit=TRUE)
aeq(p1[keep], p2)
aeq(p1, p3$fit)
aeq(p1[keep], p4$fit)
aeq(p3$se.fit[keep], p4$se.fit)

p1 <- predict(myfit, type='terms')
p2 <- predict(myfit, type="terms",newdata=lung[keep,])
p3 <- predict(myfit, type='terms', se.fit=T)
p4 <- predict(myfit, type="terms",newdata=lung[keep,], se.fit=T)
aeq(p1[keep,], p2)
aeq(p1, p3$fit)
aeq(p1[keep,], p4$fit)
aeq(p3$se.fit[keep,], p4$se.fit)

