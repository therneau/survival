cipoisson <- function(k, time=1, p=0.95, method=c('exact', 'anscombe')
{
  nn <- max(length(k), length(time), length(p))
  if (nn>1) {
      k <- rep(k, length=nn)
      time <- rep(time, length=nn)
      p <- rep(p, length=nn)
  }
  p <- (1-p)/2
  method <- match.arg(method)

  if (method=='exact') {
      dummy1 <- ifelse(k==0, 1, k)  #avoid an error message of qgamma
      lower <-  ifelse(k==0, 0, qgamma(p, dummy1))
      upper <- qgamma(1-p, k+1)
  } else if (method=='anscombe'){ # anscombe's method
      upper <- (sqrt(k + 7/8) - qnorm(p)/2)^2
      lower <- (sqrt(k - 1/8) + qnorm(p)/2)^2
  }
  else  stop("Invalid method")

  # The summary.pyears routine sometimes calls this with time=0
  if (any(time<=0)) {
      lower <- ifelse(time<=0, NA, lower)
      upper <- ifelse(time<=0, NA, upper)
  }

  if (nn==1) c(lower=lower, upper=upper)/time
  else {
      temp <- cbind(lower=lower, upper=upper)/time
      if (is.array(k)) {
         if (is.null(dd <- dimnames(k))) 
             array(temp, dim=c(dim(k), 2),
                   dimnames=list(NULL, NULL, c("lower", "upper")))
         else array(temp, dim=c(dim(k), 2),
                    dimnames=list(dd, c("lower", "upper")))
      }
      else temp
}
