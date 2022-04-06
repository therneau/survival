#
# Test bed for iteration of the model
#
itest <- function(y, x, init, offset=0, iter=20, pow=3) {
    # y = a survival object
    # x = predictor matrix, usually centered and scaled
    # init = initial values

    nvar <- ncol(x)
    beta <- matrix(0., iter, nvar)
    loglik <- dlog <- lambda <- ratio <- double(iter)
    inf <- matrix(FALSE, iter, 3)

    if (!missing(init)) {
        if (length(init) != nvar) stop("wrong length for inite")
        beta[1,] <- init
    }              
    if (missing(offset)) xoff <- rep(0., nrow(x))
    else if (length(offset) == nrow(x)) xoff <- offset
    else stop("wrong length for offset")

    ndeath <- sum(y[,ncol(y)])
    good <- 1   # last good iteration
 
    for (i in 1:iter) {
        eta <- c(x %*% beta[i,])
        if (!all(is.finite(exp(eta))) || all(exp(eta)==0)) {
            inf[i,1] <- TRUE
            loglik[i] <- loglik[i-1]
        }   
        else {
            tfit <- coxph(y ~ x + offset(xoff), iter=0, init=beta[i,], x=TRUE)
            inf[i,2] <- any(diag(tfit$var) ==0)
            inf[i,3] <- !is.finite(tfit$loglik[1])

            loglik[i] <- tfit$loglik[1]
            if (i>1) ratio[i] <- (loglik[i]- loglik[good])/ dlog[good]
            else ratio[i] <- 1
        
            dt <- coxph.detail(tfit)
            if (nvar > 1) {
                u <- colSums(dt$score)
                imat <- apply(dt$imat, 1:2, sum)
                }
             else {
                u <- sum(dt$score)
                imat <- sum(dt$imat)
             }       
        }
        if (nvar > 1) {
            imat2 <- imat + lambda[i]*diag(nvar)
            step <- drop(ginv(imat2) %*% u)
            dlog[i] <- drop(step%*% u -  step %*%imat %*% step/2)
        }   
        else {
            step <- u/(imat+lambda[i])
            dlog[i] <- step*u - 0.5* step^2*imat
            if (dlog[i] < 0) browser()
        }       
 
        if (!(any(inf[i,])) && i>1 && loglik[i] > loglik[good]) good <- i
        if (i< iter) {
            beta[i+1,] <- beta[good,] + step
                
            if (ratio[i] < .25 || any(inf[i,])) {
                if (lambda[i] >0) {lambda[i+1] <- lambda[i] * pow}
                else lambda[i+1] <- ndeath/pow
            }
            else if (ratio[i] > .75) lambda[i+1] <- lambda[i]/pow
        }
    }        
   
   data.frame(beta=beta, inf=inf, loglik=loglik, lambda=lambda, ratio=ratio, 
        dlog=dlog)
}


# Colon example
ctest <- with(colon, itest(Surv(time, status), scale(colon$time)))

# model4you, which has a horrible initial iteration
set.seed(1212)
n <- 90
d1 <- data.frame(y = abs(rnorm(n) +5) + .5, x= 1:n -10,
                    trt= rep(1:3, each=n/3))
mtest <- itest(Surv(d1$y), as.matrix(d1$trt), offset=d1$x, iter=20)

#overflow data
adata <- read.csv('overflow.csv')
atest <- with(adata, itest(Surv(start, end, status), scale(v2), iter=20))
atest2 <- with(adata, itest(Surv(start, end, status), as.matrix(v2), iter=20))
afit2 <- coxph(Surv(start, end, status) ~ v2, init=.0035, adata)
afit <- coxph(Surv(start, end, status) ~ v2, adata)

