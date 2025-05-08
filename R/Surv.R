#
# Package up surivival type data as a structure
#
Surv <- function(time, time2, event,
	      type=c('right', 'left', 'interval', 'counting', 'interval2'),
		       origin=0) {

    if (missing(time)) stop ("Must have a time argument")
    if (inherits(time ,"difftime")) time <- unclass(time)
    if (!missing(time2) && inherits(time2, "difftime")) time2 <-as.numeric(time2)
    if (!is.numeric(time)) stop ("Time variable is not numeric")
    nn <- length(time)

    # ng = number of the first 3 arguments that are present
    ng <- (!missing(time)) + (!missing(time2)) + (!missing(event)) 
    # The logic below uses "ng" throughout; why not use "missing(time2)"
    # and missing(event) instead?  Because we want to assume that 
    # "Surv(a,b)" has the variable b matched to event rather than time2.
    #
    mtype <- match.arg(c(type, 'mstate'))
    if (type== "mstate") 
        warning("type= 'mstate' is depricated, use a factor variable as status")
    
    # If type is missing or it is "mstate", I need to figure out for myself
    #  whether I have (time, time2, status) or (time, status) data
    if (missing(type) || mtype=="mstate") {
	if (ng==1 || ng==2) type <- 'right'
	else if (ng==3)     type <- 'counting'
        else stop ("No time variable!") # no time variable at all!
    }
    else {
        type <- mtype
	if (ng!=3 && (type=='interval' || type =='counting'))
		stop("Wrong number of args for this type of survival data")
	if (ng!=2 && (type=='right' || type=='left' ||  type=='interval2'))
		stop("Wrong number of args for this type of survival data")
	}

    if (ng==1) { # only a time variable given
        if (!is.numeric(time)) stop("Time variable is not numeric")
	ss <- cbind(time=time-origin, status=1)
        type <- "right"
	}
    else if (type=='right' || type=='left') {
        if (!is.numeric(time)) stop("Time variable is not numeric")
	if (missing(event))   {
            event <- time2  # treat time2 as event
            time2 <- NULL   # force any inputAttributes to attach to "event"
        }
        if (length(event) != nn) stop ("Time and status are different lengths")
        if (mtype=="mstate" || (is.factor(event))) {
            mstat <- as.factor(event)
            status <- as.numeric(mstat) -1
            type <- "mright"
        }
	else {
            if (is.logical(event)) status <- as.numeric(event)
            else  if (is.numeric(event)) {
                who2 <- !is.na(event)
                if (max(event[who2]) ==2) status <- event -1
                else status <- event
                temp <- (status==0 | status==1)
                status <- ifelse(temp, status, NA)
                if (!all(temp[who2], na.rm=TRUE))
		    warning("Invalid status value, converted to NA")
	    }
            else stop("Invalid status value, must be logical or numeric")
        }
	ss <- cbind(time=time-origin, status=status)
    }
    else  if (type=='counting') {
	if (length(time2) !=nn) stop ("Start and stop are different lengths")
	if (length(event)!=nn) stop ("Start and event are different lengths")
        if (!is.numeric(time))  stop("Start time is not numeric")
	if (!is.numeric(time2)) stop("Stop time is not numeric")
	temp <- (time >= time2)
	if (any(temp & !is.na(temp))) {
	    time[temp] <- NA
	    warning("Stop time must be > start time, NA created")
	    }
        if (mtype=="mstate" || is.factor(event)) {
            mstat <- as.factor(event)
            status <- as.numeric(mstat) -1
            type <- "mcounting"
        }
	else {
            if (is.logical(event)) status <- as.numeric(event)
	    else  if (is.numeric(event)) {
		who2 <- !is.na(event)
		if (max(event[who2])==2) status <- event - 1
		else status <- event
		temp <- (status==0 | status==1)
		status <- ifelse(temp, status, NA)
		if (!all(temp[who2], na.rm=TRUE))
		    warning("Invalid status value, converted to NA")
		}
	    else stop("Invalid status value")
        }
	ss <- cbind(start=time-origin, stop=time2-origin, status=status)
    }

    else {  #interval censored data
	if (type=='interval2') {
	    # convert to "interval" type, infer the event code
	    if (!is.numeric(time2)) stop("Time2 must be numeric")
	    if (length(time2) !=nn) 
		    stop ("time and time2 are different lengths")
            backwards <- (!is.na(time) & !is.na(time2) & time > time2)
            # allow for infinite times (important to do the backwards check
            #  first)
            time  <- ifelse(is.finite(time), time, NA)
            time2 <- ifelse(is.finite(time2), time2, NA)
            unknown <-  (is.na(time) & is.na(time2))
 	    status <- ifelse(is.na(time),  2,
		      ifelse(is.na(time2), 0,
		      ifelse(time==time2, 1,3)))
	    time <- ifelse(status!=2, time, time2)

            if (any(backwards)) {
                warning("Invalid interval: start > stop, NA created")
                status[backwards] <- NA
            }
            if (any(unknown)) status[unknown] <- NA
	    type <- 'interval'
	    }
	else {  #check legality of event code
	    if (length(event)!=nn) 
		    stop("Time and status are different lengths")
            if (!is.numeric(event)) 
		   stop("Invalid status value, must be logical or numeric")
            temp <- (event==0 | event==1| event==2 | event==3)
            status <- ifelse(temp, event, NA)
            if (!all(temp, na.rm=TRUE))
                warning("Status must be 0, 1, 2 or 3; converted to NA")

	    if (any(event==3, na.rm=T)) {
		if (!is.numeric(time2)) stop("Time2 must be numeric")
		if (length(time2) !=nn) 
		    stop ("time and time2 are different lengths")
                temp <- (status==3 & time>time2)
                if (any(temp & !is.na(temp))) {
                    status[temp] <- NA
                    warning("Invalid interval: start > stop, NA created")
                }
            }
	    else time2 <- 1  #dummy value, time2 is never used
        }

	ss <- cbind(time1=time-origin, 
		    time2=ifelse(!is.na(status) & status==3, time2-origin, 1),
		    status=status)
    }

    # Retain any attributes of the input arguments. Originally requested
    #  by the rms package
    inputAttributes <- list()
    if (!is.null(attributes(time)))
        inputAttributes$time  <-attributes(time)
    if (!missing(time2) && !is.null(attributes(time2)))
        inputAttributes$time2 <- attributes(time2)
    if (!missing(event) && !is.null(attributes(event)))
        inputAttributes$event <- attributes(event)

    # In rare cases there are no column names, and I have discovered that
    #  people depend on them.
    cname <- dimnames(ss)[[2]]
    if (length(cname) ==0) {
        if (ncol(ss)==2) cname <- c("time", "status")
        else if (type=="counting") cname <- c("start", "stop", "status")
        else cname <- c("time1", "time2", "status")
    }
    dimnames(ss) <- list(NULL, cname)  #kill extraneous row names
                                           
    attr(ss, "type")  <- type
    if (type=="mright" || type=="mcounting") {
        states <- levels(mstat)[-1]
        if (any(is.na(states) | states=='') )
            stop("each state must have a non-blank name")
        attr(ss, "states") <- states
    }
    if (length(inputAttributes) > 0) 
        attr(ss, "inputAttributes") <- inputAttributes
    class(ss) <- 'Surv'
    ss
    }

print.Surv <- function(x, quote=FALSE, ...) {
    invisible(print(as.character.Surv(x), quote=quote, ...))
    }

as.character.Surv <- function(x, ...) {
    new <- switch(attr(x, "type"),
           "right"={
               temp <- x[,2]
               temp <- ifelse(is.na(temp), "?", ifelse(temp==0, "+",""))
               paste0(format(x[,1]), temp)
           },
           "counting"= {
               temp <- x[,3]
               temp <- ifelse(is.na(temp), "?", ifelse(temp==0, "+",""))
               paste0('(', format(x[,1]), ',', format(x[,2]), temp,
                     ']')
           },
           "left" ={
               temp <- x[,2]
               temp <- ifelse(is.na(temp), "?", ifelse(temp==0, "-",""))
               paste0(format(x[,1]), temp)
           },
           "interval"= {
               stat <- x[,3]
               temp <- c("+", "", "-", "]")[stat+1]
               temp2 <- ifelse(stat==3,
			 paste("[", format(x[,1]), ", ",format(x[,2]), sep=''),
			 format(x[,1]))
               ifelse(is.na(stat), "NA", paste0(temp2, temp))
           },
           "mright" = {  #multi-state
               temp <- x[,2]
               end <- c("+", paste(":", attr(x, "states"), sep='')) #endpoint
               temp <- ifelse(is.na(temp), "?", end[temp+1])
               paste0(format(x[,1]), temp)
           },
           "mcounting"= {
               temp <- x[,3]
               end <- c("+", paste(":", attr(x, "states"), sep='')) #endpoint
               temp <- ifelse(is.na(temp), "?", end[temp+1])
               paste0('(', format(x[,1]), ',', format(x[,2]), temp,
                     ']')
           })
    names(new) <- rownames(x)
    new
}


"[.Surv" <- function(x, i, j, drop=FALSE) {
    # If only 1 subscript is given, the result will still be a Surv object,
    #   and the drop argument is ignored.
    # I would argue that x[3:4,,drop=FALSE] should return a matrix, since
    #  the user has implicitly specified that they want a matrix.
    #  However, [.dataframe calls [.Surv with the extra comma; its
    #  behavior drives the choice of default.
    if (missing(j)) {
        xattr <- attributes(x)
        x <- unclass(x)[i,, drop=FALSE] # treat it as a matrix: handles dimnames
        attr(x, 'type') <- xattr$type
        if (!is.null(xattr$states)) attr(x, "states") <- xattr$states
        if (!is.null(xattr$inputAttributes)) {
            # If I see "names" subscript it, leave all else alone
            attr(x, 'inputAttributes') <- 
                lapply(xattr$inputAttributes, function(z) {
                       if (any(names(z)=="names")) z$names <- z$names[i]
                       z
                   })
        }
        class(x) <- "Surv"  #restore the class
        x
    }
    else { # return  a matrix or vector
	class(x) <- 'matrix'
       	NextMethod("[")
    }
}

is.na.Surv <- function(x) {
    as.vector(rowSums(is.na(unclass(x))) >0)
    }

Math.Surv <- function(...)  stop("Invalid operation on a survival time")
Ops.Surv  <- function(...)  stop("Invalid operation on a survival time")
Summary.Surv<-function(...) stop("Invalid operation on a survival time")

# The Ops.Surv method could in theory define == and >, to allow sorting
#  but I've left them out since it is the xtfrm method that explicitly
#  is used for this.  For (start, stop) data we order by event within
#  ending time.  Start time is included as a last index, but it is not
#  clear that we need to do so.
xtfrm.Surv <- function(x) {
    if (attr(x, 'type') == "interval") {
        temp <- ifelse(x[,3]==3, (x[,1] + x[,2])/2, x[,1])
        index <- order(temp, match(x[,3], c(2,1,3,0)))
        }
    else if (attr(x, 'type')== "left") index <- order(x[,1], x[,2])
    else if (ncol(x)==2) index <- order(x[,1], x[,2]==0, x[,2]) # censor last
    else index <- order(x[,2], x[,3]==0, x[,3], x[,1]) # ending time, stat, start
    temp <- integer(nrow(x))
    temp[index] <- seq.int(nrow(x))
    temp[is.na(x)] <- NA
    temp
}

is.Surv <- function(x) inherits(x, 'Surv')
as.matrix.Surv <- function(x, ...) {
    y <- unclass(x)
    attr(y, "type") <- NULL
    attr(y, "states") <- NULL
    attr(y, "inputAttributes") <- NULL
    y
}

# You can't do length without names
# and names doesn't pay attention to my definition of length:
# we need to map to rownames instead

length.Surv <- function(x) nrow(x)
"names<-.Surv" <- function(x, value) {
    rownames(x) <- value
    x
}
names.Surv <- function(x) rownames(x)

format.Surv <- function(x, ...) format(as.character.Surv(x), ...)
as.data.frame.Surv <- function(x, ...) as.data.frame.model.matrix(x, ...)

# all sorts of methods for Surv, caused by searching for every case of
#  UseMethod in the standard libraries

# package:utils methods
tail.Surv <- function(x, ...) 
    x[tail(1:nrow(x), ...)]     

head.Surv <- function(x, ...)
    x[head(1:nrow(x), ...)]

# packge:graphics.  All try to give a nicer failure message
plot.Surv <- function(x, ...)
    plot(survfit(x ~1), ...)

barplot.Surv <- function(height, ...)
    stop("not defined for a Surv object")
density.Surv <- function(x, ...)
    stop("method not defined for a Surv object")
hist.Surv <- function(x, ...)
    stop("method not defined for a Surv object")
identify.Surv <- function(x, ...)
    stop("method not defined for a Surv object")
image.Surv <- function(x, ...)
    stop("method not defined for a Surv object")
lines.Surv <- function(x, ...)
    stop("method not defined for a Surv object")
points.Surv <- function(x, ...)
    stop("method not defined for a Surv object")
pairs.Surv <- function(x, ...)
    stop("method not defined for a Surv object")
text.Surv <- function(x, ...)
    stop("method not defined for a Surv object")

# package:base methods
anyDuplicated.Surv <- function(x, ...) anyDuplicated(as.matrix(x), ...)
duplicated.Surv    <- function(x, ...) duplicated(as.matrix(x), ...)
rev.Surv <- function(x) x[rev(1:nrow(x))]
unique.Surv  <- function(x, ...)
    x[!duplicated(as.matrix(x), ...)]

c.Surv <- function(...) {
    slist <- list(...)
    if (!all(sapply(slist, function(x) inherits(x, "Surv")))) 
        stop("all elements must be of class Surv")
    types <- sapply(slist, function(x) attr(x, "type"))
    if (!all(types == types[1]))
        stop("all elements must be of the same Surv type")

    if (types[[1]] %in% c("mright", "mcounting")) {
        states <- lapply(slist, function(x) attr(x, 'states'))
        if (any(diff(sapply(states, length))!=0))
            stop("all elements must have the same list of states")
        if (!all(sapply(states, function(x) all.equal(x, states[[1]]))))
            stop("all elements must have the same list of states")
        }
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
#rbind.Surv <- function(...) {
#    dotlist <- list(...)
#    if (all(sapply(dotlist, is.Surv))) do.call("c.Surv", dotlist)
#    else do.call("rbind", lapply(dotlist, function(x)
#        if (is.Surv(x)) as.matrix(x) else x))
#    }
#
#cbind.Surv <- function(...) 
#    do.call("cbind", lapply(list(...), 
#                        function(x) if (is.Surv(x)) as.matrix(x) else x))
#}

rep.Surv <- function(x, ...) {
    index <- rep(1:nrow(x), ...)
    x[index]
    }
rep.int.Surv <- function(x, ...) {
    index <- rep.int(1:nrow(x), ...)
    x[index]
    }
rep_len.Surv <- function(x, ...) {
    index <- rep_len(1:nrow(x), ...)
    x[index]
    }
t.Surv <- function(x) t(as.matrix(x))

as.logical.Surv <- function(x, ...)
    stop("invalid operation on a survival time")

# removed 2024-06-02, make Surv act like a matrix for this op
#as.integer.Surv <- function(x, ...) {
#    nc <- ncol(x)
#    x[,-nc] <- as.integer(x[,-nc])
#    if (nc==3 && any(x[,1] >= x[,2]))
#        stop("invalid survival time created")
#    x
#}

# per the help file for as.numeric, this should have been as.double
#  so never worked anyway
#as.numeric.Surv <- function(x, ...) {
#    nc <- ncol(x)
#    x[,-nc] <- as.numeric(x[, -nc])
#    x
#}

mean.Surv <-function(x, ...)
    stop("a mean method has not been defined for Surv objects")

median.Surv <- function(x, na.rm=FALSE, ...)
    quantile(x, probs=0.5, na.rm= na.rm, ...)

quantile.Surv <- function(x, probs, na.rm=FALSE, ...) {
    if (!na.rm && any(is.na(x)))
        stop("missing values and NaN's not allowed if 'na.rm' is FALSE")
    if (attr(x, "type") %in% c("mright", "mcounting"))
        stop("quantile method not defined for multiple-endpoint Surv objects")
    fit <- survfit(x~1)
    quantile.survfit(fit, probs, ...)
}   

# these make sense but aren't S3 methods
# sd, IQR, mad, cov, cor

levels.Surv <- function(x) attr(x, "states")
