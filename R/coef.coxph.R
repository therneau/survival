coef.coxphms <- function(object, type=c("vector", "matrix"), ...) {
    type <- match.arg(type)
    if (type=="matrix") {
        cmap2 <- object$cmap[-1,, drop=FALSE]
        cmat <- 0*cmap2  # all the right names
        cmat[cmap2>0] <- object$coefficient[cmap2]
        attr(cmat, "states") <- object$states
        cmat
    }
    else NextMethod(object, ...)
}

        
