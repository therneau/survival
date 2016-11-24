# License: 
# 
# Copyright 2004 Mayo Foundation for Medical Education and Research. 
# 
# This program is free software; you can redistribute it and/or modify it under the terms of 
# the GNU General Public License as published by the Free Software Foundation; either 
# version 2 of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT ANY 
# WARRANTY; without even the implied warranty of MERCHANTABILITY or 
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
# more details.


##
## $Id: cipoisson.q,v 1.1.1.1 2011/02/24 14:16:27 sinnwell Exp $
#
# Exact Poisson confidence limits
#   Based on equation 10.10 of Feller volume 1, chapter 6
# The other approximations make use of the continuity correction

#' Confidence limits for the Poisson
#' 
#' Give exact confidence limits for a Poisson parameter.
#' 
#' @param k Number of successes
#' @param time Total time on trial
#' @param p Probability level for the (two-sided) interval
#' @param method The method of approximation. The approximate method is supplied for comparison purposes.
#' @return Either a vector containing the lower and upper limits, or a two column matrix of lower and upper limits.
#' @details The exact method is based on equation 10.10 of Feller, which relates poisson probabilities to tail area of the gamma distribution.
#'   The Anscombe approximation is based on the fact that sqrt(k + 3/8) is has a nearly constant variance of 1/4,
#'   along with a continuity correction.
#' 
#' @references W.F. Feller (1950). An Introduction to Probability Theory and its Applications, Volume 1, Chapter 6, Wiley.
#' 
#' F.J. Anscombe (1949). Transformations of Poisson, binomial and negative-binomial data. Biometrika, vol 35, p246-254.
#' 
#' @examples
#' cipoisson(4) # 95\% confidence limit 
#' # lower    upper  
#' # 1.089865 10.24153 
#' ppois(4, 10.24153)     #chance of seeing 4 or fewer events with large rate  
#' # [1] 0.02500096 
#' 1-ppois(3, 1.08986)    #chance of seeing 4 or more, with a small rate 
#' # [1] 0.02499961
#' 
#' @seealso \code{\link[stats]{ppois}}, \code{\link[stats]{qpois}}
#' @export

cipoisson <- function(k, time=1, p=0.95, method=c('exact', 'anscombe'))
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
  } else stop("Invalid method")

  if (nn==1) c(lower=lower, upper=upper)/time
  else       cbind(lower=lower, upper=upper)/time
}
