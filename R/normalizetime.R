#
# Create time values such that tiny differences are treated as a tie
#  The actions and arguments are the same as all.equal
#
normalizetime <- function(x, tolerance = sqrt(.Machine$double.eps)) {
    if (is.Surv(x)) y <- sort(x[, -ncol(x)])
    else y <- sort(x)

    dy <- diff(y)
    tied <- ((dy <=tolerance) |( (dy/ mean(abs(y)) <=tolerance)))
    newy <- y[c(!tied, TRUE)]
    
    if (is.Surv(x)) {
        x[, -ncol(x)] <- findInterval(x[, -ncol(x)], newy)
        attr(x, 'utime') <- newy
        x
    } else {
        z <- findInterval(x, newy)
        attributes(z) <- attributes(x)
        attr(z, 'utime') <- newy
        class(z) <- "normalizetime"
        z}
}
