library(survival)
library(splines)

# the nsk function should give the same solution as ns, but with a different
#  parameterization
#
xx <- runif(500, 1, 100)
yy <- 10*log(xx) + rnorm(500, 0, 2)
tdata <- data.frame(xx=xx, yy=yy)
fit1 <- lm(yy ~ ns(xx, df=4), tdata, model=TRUE)
fit2 <- lm(yy ~ nsk(xx, df=4, b=0), tdata)
all.equal(predict(fit1), predict(fit2))  # same solution

xattr <- attributes(fit1$model[[2]])
allknots <- sort(c(xattr$knots, xattr$Boundary.knots)) # knots that were used
pred.knot <- predict(fit1, newdata=list(xx=allknots))
all.equal(pred.knot[-1] - pred.knot[1], coef(fit2)[-1],
          check.attributes = FALSE)

