# The density, quantile, and CDF functions for those distributions
#   supported by survreg
#
dsurvreg <- function(x, mean, scale=1, distribution='weibull', parms) {
    dist <- survreg.distributions[[casefold(distribution)]]
    if (is.null(dist)) stop("Distribution not found")

    if (!is.null(dist$trans)) {
	dx <- dist$dtrans(x)
	x  <- dist$trans(x)
	x <- (x-mean)/scale
	dist <- survreg.distributions[[dist$dist]]
	y <- dist$density(x, parms)[,3]
	y *dx / scale
	}
    else {
	x <- (x-mean)/scale
	y <- dist$density(x, parms)[,3]
	y/ scale
	}
    }

psurvreg <- function(q, mean, scale=1, distribution='weibull', parms) {
    dist <- survreg.distributions[[casefold(distribution)]]
    if (is.null(dist)) stop("Distribution not found")

    if (!is.null(dist$trans)) {
	q  <- dist$trans(q)
	q <- (q-mean)/scale
	dist <- survreg.distributions[[dist$dist]]
	dist$density(q, parms)[,1]
	}
    else {
	q <- (q-mean)/scale
	dist$density(q, parms)[,1]
	}
    }

qsurvreg <- function(p, mean, scale=1, distribution='weibull', parms) {
    dist <- survreg.distributions[[casefold(distribution)]]
    if (is.null(dist)) stop("Distribution not found")

    if (!is.null(dist$trans)) {
	d2 <- survreg.distributions[[dist$dist]]
	x  <- d2$quantile(p, parms)
	dist$itrans(x*scale + mean)
	}
    else {
	x <- dist$quantile(p, parms)
	x*scale + mean
	}
    }

rsurvreg <- function(n, mean, scale=1, distribution='weibull', parms) {
    if (missing(parms)) 
        qsurvreg(runif(n), mean, scale, distribution)
    else
        qsurvreg(runif(n), mean, scale, distribution, parms) 
}
