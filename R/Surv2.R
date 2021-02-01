#
# Package up surivival type data as a structure
#  This variant is for time-course data
#
Surv2 <- function(time, event, repeated=FALSE) {
    if (missing(time)) stop ("must have a time argument")
    if (inherits(time ,"difftime")) time <- unclass(time)
    if (!is.numeric(time)) stop ("Time variable is not numeric")
    nn <- length(time)
    if (!is.logical(repeated) || length(repeated) != 1)
        stop("invalid value for repeated option")

    if (missing(event)) stop("must have an event argument")
    if (length(event) != nn) stop ("Time and event are different lengths")
    event <- as.factor(event)
    states <- levels(event)[-1]
    status <- as.integer(event) -1L # usually time is not integer, but
    ss <- cbind(time=time, status=status) # sometimes it is
 
    # In rare cases there are no column names, and I have discovered that
    #  people depend on them.
    cname <- dimnames(ss)[[2]]
    if (length(cname) ==0) cname <- c("time", "status")
    dimnames(ss) <- list(NULL, cname)  #kill extraneous row names
                                           
    if (any(is.na(states) | states=='') )
        stop("each state must have a non-blank name")
    attr(ss, "states") <- states
    attr(ss, "repeated") <- repeated
    class(ss) <- 'Surv2'
    ss
    }

print.Surv2 <- function(x, quote=FALSE, ...) {
    invisible(print(as.character.Surv2(x), quote=quote, ...))
    }

as.character.Surv2 <- function(x, ...) {
    states <- attr(x, "states")
    if (is.null(states)) {
        temp <- x[,2]
        temp <- ifelse(is.na(temp), "?", ifelse(temp==0, "+",""))
        new <- paste0(format(x[,1]), temp)
    } else {
        temp <- x[,2]
        end <- c("+", paste(":", states, sep='')) #endpoint
        temp <- ifelse(is.na(temp), "?", end[temp+1])
        new <- paste0(format(x[,1]), temp)
    }
    names(new) <- rownames(x)
    new
}


"[.Surv2" <- function(x, i, j, drop=FALSE) {
    # If only 1 subscript is given, the result will still be a Surv2 object,
    #   and the drop argument is ignored.
    # I would argue that x[3:4,,drop=FALSE] should return a matrix, since
    #  the user has implicitly specified that they want a matrix.
    #  However, [.dataframe calls [.Surv with the extra comma; its
    #  behavior drives the choice of default.
    if (missing(j)) {
        xattr <- attributes(x)
        x <- unclass(x)[i,, drop=FALSE] # treat it as a matrix: handles dimnames
        if (!is.null(xattr$states)) attr(x, "states") <- xattr$states
        attr(x, "repeated") <- xattr$repeated
        class(x) <- "Surv2"  #restore the class
        x
    }
    else { # return  a matrix or vector
	class(x) <- 'matrix'
       	NextMethod("[")
    }
}

is.na.Surv2 <- function(x) {
    as.vector(rowSums(is.na(unclass(x))) >0)
    }

Math.Surv2 <- function(...)  stop("Invalid operation on a survival time")
Ops.Surv2  <- function(...)  stop("Invalid operation on a survival time")
Summary.Surv2<-function(...) stop("Invalid operation on a survival time")

# The Ops.Surv method could in theory define == and >, to allow sorting
#  but I've left them out since it is the xtfrm method that explicitly
#  is used for this.  For (start, stop) data we order by event within
#  ending time.  Start time is included as a last index, but it is not
#  clear that we need to do so.
xtfrm.Surv2 <- function(x) {
    index <- order(x[,1], x[,2]==0, x[,2]) # censor last
    temp <- integer(nrow(x))
    temp[index] <- seq.int(nrow(x))
    temp[is.na(x)] <- NA
    temp
}

is.Surv2 <- function(x) inherits(x, 'Surv2')
as.matrix.Surv2 <- function(x, ...) {
    y <- unclass(x)
    attr(y, "states") <- NULL
    attr(y, "repeated") <- NULL
    y
}

# You can't do length without names
# and names doesn't pay attention to my definition of length:
# we need to map to rownames instead

length.Surv2 <- function(x) nrow(x)
"names<-.Surv2" <- function(x, value) {
    rownames(x) <- value
    x
}
names.Surv2 <- function(x) rownames(x)

format.Surv2 <- function(x, ...) format(as.character.Surv2(x), ...)
as.data.frame.Surv2 <- function(x, ...) as.data.frame.model.matrix(x, ...)

# all sorts of methods for Surv, caused by searching for every case of
#  UseMethod in the standard libraries

# package:utils methods
tail.Surv2 <- function(x, ...) 
    x[tail(1:nrow(x), ...)]     

head.Surv2 <- function(x, ...)
    x[head(1:nrow(x), ...)]

# packge:graphics.  All try to give a nicer failure message
plot.Surv2 <- function(x, ...)
    stop("method not defined for a Surv2 object")

barplot.Surv2 <- function(height, ...)
    stop("not defined for a Surv2 object")
density.Surv2 <- function(x, ...)
    stop("method not defined for a Surv2 object")
hist.Surv2 <- function(x, ...)
    stop("method not defined for a Surv2 object")
identify.Surv2 <- function(x, ...)
    stop("method not defined for a Surv2 object")
image.Surv2 <- function(x, ...)
    stop("method not defined for a Surv2 object")
lines.Surv2 <- function(x, ...)
    stop("method not defined for a Surv2 object")
points.Surv2 <- function(x, ...)
    stop("method not defined for a Surv2 object")
pairs.Surv2 <- function(x, ...)
    stop("method not defined for a Surv2 object")
text.Surv2 <- function(x, ...)
    stop("method not defined for a Surv2 object")

# package:base methods
anyDuplicated.Surv2 <- function(x, ...) anyDuplicated(as.matrix(x), ...)
duplicated.Surv2   <- function(x, ...) duplicated(as.matrix(x), ...)
rev.Surv2 <- function(x) x[rev(1:nrow(x))]
unique.Surv2 <- function(x, ...)
    x[!duplicated(as.matrix(x), ...)]

c.Surv2 <- function(...) {
    slist <- list(...)
    if (!all(sapply(slist, function(x) inherits(x, "Surv2")))) 
        stop("all elements must be of class Surv2")

    states <- lapply(slist, function(x) attr(x, 'states'))
    if (any(diff(sapply(states, length))!=0))
        stop("all elements must have the same list of states")
    if (!is.null(states[[1]]) && 
                 !all(sapply(states, function(x) all.equal(x, states[[1]]))))
            stop("all elements must have the same list of states")
  
    new <- do.call("rbind", lapply(slist, as.matrix))
    att1 <- attributes(slist[[1]])
    att1 <- att1[is.na(match(names(att1), c("dim","dimnames")))]
    attributes(new) <- c(attributes(new)[c('dim', 'dimnames')], att1)
    new
    }


# The cbind/rbind methods cause more trouble than they solve
# The problem is when one is called with mixed arguments, e.g.
#      cbind(Surv(1:4), data.frame(x=6:9, z=c('a', 'b', 'a', 'a'))
# R will call cbind.Surv for cbind(Surv(1:4), Surv(2:5)), giving a matrix.
# In the above, however, cbind.Surv is never called, but the \emph{presence}
#    of a cbind.Surv method messes up the default behavior, see the
#    'Dispatch' section of help('cbind').  The result becomes a matrix of
#    lists rather than a dataframe.
#

rep.Surv2 <- function(x, ...) {
    index <- rep(1:nrow(x), ...)
    x[index]
    }
rep.int.Surv2 <- function(x, ...) {
    index <- rep.int(1:nrow(x), ...)
    x[index]
    }
rep_len.Surv2 <- function(x, ...) {
    index <- rep_len(1:nrow(x), ...)
    x[index]
    }
t.Surv2 <- function(x) t(as.matrix(x))

as.logical.Surv2 <- function(x, ...)
    stop("invalid operation on a survival time")
as.integer.Surv2 <- function(x, ...) {
    x[,1] <- as.integer(x[,1])
    x
}
as.numeric.Surv2 <- function(x, ...) {
    x[,1] <- as.numeric(x[,1])
    x
}

mean.Surv2 <-function(x, ...)
    stop("a mean method has not been defined for Surv2 objects")

median.Surv2 <- function(x, ...)
    stop("median method has not been defined for Surv2 objects")

quantile.Surv2 <- function(x, probs, na.rm=FALSE, ...) {
    if (!na.rm && any(is.na(x)))
        stop("missing values and NaN's not allowed if 'na.rm' is FALSE")
    stop("quantile method not defined for Surv2 objects")
}   

# these make sense but aren't S3 methods
# sd, IQR, mad, cov, cor

levels.Surv2 <- function(x) attr(x, "states")
