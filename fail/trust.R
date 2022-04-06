# A look at trust region ideas
#  First the quadratic approx underlying Newton Raphson
aeq <- function(x, y, ...)  all.equal(as.vector(x), as.vector(y), ...)

quad <- function(first, second) {
    delta <- solve(second, first)
    yhat <-  sum(first*delta) - delta %*% second %*% delta/2
    list(delta=delta, yhat=c(yhat))
}


test0 <- coxph(Surv(time, status) ~ ph.ecog + age + sex, lung, iter=0)
test1 <- coxph(Surv(time, status) ~ ph.ecog + age + sex, lung, iter=1)
dt <- coxph.detail(test0)

temp <- quad(colSums(dt$score), apply(dt$imat, 1:2, sum))
aeq(temp$delta, coef(test1) - coef(test0))
temp$yhat/diff(test1$loglik)  # accuracy of the approximation
         
# Look at the first iterations of an infinite failure
c0 <- coxph(Surv(time, status) ~ I(-time/1000), colon, iter=0)
dt <- coxph.detail(c0)

beta <- seq(0, 4, length=50)
ll <- matrix(0, 50, 2)
for (i in 1:50) {
    ll[i,1] <- coxph(Surv(time, status) ~ I(-time/1000), colon, iter=0, 
               init=beta[i])$loglik[1] - c0$loglik[1]
    ll[i,2] <- beta[i]*sum(dt$score) - beta[i]^2*sum(dt$imat)/2
}
matplot(beta, ll, type='l')

c1 <- coxph(Surv(time, status) ~ I(-time/1000), colon, iter=1)
dt <- coxph.detail(c1)

beta <- seq(2, 6, length=50)
ll <- matrix(0, 50, 2)
for (i in 1:50) {
    ll[i,1] <- coxph(Surv(time, status) ~ I(-time/1000), colon, iter=0, 
               init=beta[i])$loglik[1] - c1$loglik[2]
    d <- beta[i] - c1$coef
    ll[i,2] <- d*sum(dt$score) - d^2*sum(dt$imat)/2
}
matplot(beta, ll, type='l')

cmat <- matrix(0, 9, 3)
for (i in 0:8) {
    cfit <- coxph(Surv(time, status) ~ time, colon, iter=i)
    cdt <- coxph.detail(cfit)
    qtemp <- quad(sum(cdt$score), sum(cdt$imat))
    cmat[i+1, 1] <- qtemp$yhat
    cmat[i+1, 2] <- cfit$coef
    cmat[i+1, 3] <- cfit$loglik[2]
}

matplot(1:8, cbind(cmat[1:8, 1], diff(cmat[,3])), 
        xlab='Iter', ylab="Change")

# The model4you error
set.seed(1212)
n <- 90
d1 <- data.frame(y = abs(rnorm(n) +5) + .5, x= 1:n -10,
                    trt= rep(1:3, each=n/3))

m0 <- coxph(Surv(y) ~ trt + offset(x), d1, iter=0)
dm <- coxph.detail(m0)
qm <- quad(sum(dm$score), sum(dm$imat))

qm2 <- quad(sum(dm$score) -1, sum(dm$imat) +1)
test2 <- coxph(Surv(y) ~ trt + offset(x), d1, iter=0, init=qm2$delta)
qm3 <- quad(sum(dm$score) -2, sum(dm$imat) +2)
qm4 <- quad(sum(dm$score) -4, sum(dm$imat) +4)
