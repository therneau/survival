#
# Create time values such that tiny differences are treated as a tie
#  The actions and tolerance are the same as all.equal
#
normalizetime <- function(x, replace=TRUE,
                          tolerance = sqrt(.Machine$double.eps)) {
    if (is.Surv(x)) y <- sort(unique(c(x[, -ncol(x)])))
    else y <- sort(unique(x))
    y <- y[is.finite(y)]  #someone may hand us an INF

    dy <- diff(y)
    tied <- ((dy <=tolerance) |( (dy/ mean(abs(y)) <=tolerance)))
    if (!any(tied)) return(x)   # all values are unique

    cuts <- y[c(TRUE, !tied)]
    if (is.Surv(x)) {
        z <- findInterval(x[, -ncol(x)], cuts)
        if (replace) {
            z <- matrix(c(cuts[z], as.integer(x[,ncol(x)])), ncol=ncol(x))
            attributes(z) <- attributes(x)
        }
        else {
            z <- matrix(c(z, as.integer(x[,ncol(x)])), ncol=ncol(x))
            attributes(z) <- attributes(x)
            attr(z, 'utime') <-  unname(cuts)
        }
    } else {
        z <- findInterval(x, cuts)
        if (replace) {
            z <- cuts[z]
            attributes(z) <- attributes(x)
        }
        else {
            attributes(z) <- attributes(x)
            attr(z, 'utime') <- unname(cuts)
        }
    }
    z
}
