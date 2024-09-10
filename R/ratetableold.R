# These functions are depricated
#  Only relsurv uses them, and I'm working on that
ratetable <- function(...) {
    datecheck <- function(x) 
        inherits(x, c("Date", "POSIXt", "date", "chron"))

    args <- list(...)
    nargs <- length(args)
    ll <- sapply(args, length)
    n <- max(ll)  # We assume this is the dimension of the user's data frame
    levlist <- vector("list", nargs)
    x <- matrix(0,n,nargs)
    dimnames(x) <- list(1:n, names(args))
    isDate <- sapply(args, datecheck)

    for (i in 1:nargs) {
        if (ll[i] ==1) args[[i]] <- rep(args[[i]], n)
        else if (ll[i] != n) 
            stop(gettextf("Aguments do not all have the same length (arg %d)", i))

	# In Splus cut and tcut produce class 'category'
        if (inherits(args[[i]], 'cateogory') || is.character(args[[i]]))
                args[[i]] <- as.factor(args[[i]])
        if (is.factor(args[[i]])) {
            levlist[[i]] <- levels(args[[i]])
            x[,i] <- as.numeric(args[[i]]) # the vector of levels
            }
        else x[,i] <- ratetableDate(args[[i]]) 
 	}
    attr(x, "isDate") <- isDate
    attr(x, "levlist")   <- levlist
    class(x) <- 'ratetable2'
    x
    }

# The two functions below should only be called internally, when missing
#   values cause model.frame to drop some rows
is.na.ratetable2 <- function(x) {
    attributes(x) <- list(dim=dim(x))
    as.vector((1 * is.na(x)) %*% rep(1, ncol(x)) >0)
    }
"[.ratetable2" <- function(x, rows, cols, drop=FALSE) {
    if (!missing(cols)) {
       stop("This should never be called!")
       }
    aa <- attributes(x)
    attributes(x) <- aa[c("dim", "dimnames")]
    y <- x[rows,,drop=FALSE]
    attr(y,'isDate') <- aa$isDate
    attr(y,'levlist')   <- aa$levlist
    class(y) <- 'ratetable2'
    y
    }
