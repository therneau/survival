#
# Package up surivival type data as a structure
#
Surv <- function(time, time2, event,
	      type=c('right', 'left', 'interval', 'counting', 'interval2',
                     "mstate"),
		       origin=0) {

    if (missing(time)) stop ("Must have a time argument")
    if (!is.numeric(time)) stop ("Time variable is not numeric")
    nn <- length(time)

    # ng = number of the first 3 arguments that is present
    ng <- (!missing(time)) + (!missing(time2)) + (!missing(event)) 
    # The logic below uses "ng" throughout; why not use "missing(time2)"
    # and missing(event) instead?  Because we want to assume that 
    # "Surv(a,b)" has the variable b matched to event rather than time2.
    #
    mtype <- match.arg(type)

    # If type is missing or it is "mstate", I need to figure out for myself
    #  whether I have (time1, time2, status) or (time, status) data
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
	if (missing(event))   event <- time2  # treat time2 as event
        if (length(event) != nn) stop ("Time and status are different lengths")
        if (mtype=="mstate" || (is.factor(event) && length(levels(event))>2)) {
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
        if (mtype=="mstate" || (is.factor(event) && length(levels(event))>2)) {
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
		    stop ("Time1 and time2 are different lengths")
	    status <- ifelse(is.na(time), 2,
		      ifelse(is.na(time2),0,
		      ifelse(time==time2, 1,3)))
	    time <- ifelse(status!=2, time, time2)
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
		    stop ("Time1 and time2 are different lengths")
		}
	    else time2 <- 1  #dummy value, time2 is never used
	    }

	temp <- (status==3 & time>time2)
	if (any(temp & !is.na(temp))) {
	    time[temp] <- NA
	    warning("Invalid interval: start > stop, NA created")
	    }

	ss <- cbind(time1=time-origin, 
		    time2=ifelse(!is.na(status) & status==3, time2-origin, 1),
		    status=status)
	}

    dimnames(ss) <- list(NULL, dimnames(ss)[[2]]) #kill any tag-along row names
    attr(ss, "type")  <- type
    if (type=="mright" || type=="mcounting") attr(ss, "states") <- levels(mstat)[-1]
    class(ss) <- 'Surv'
    ss
    }

print.Surv <- function(x, quote=FALSE, ...) {
    invisible(print(as.character.Surv(x), quote=quote, ...))
    }

as.character.Surv <- function(x, ...) {
    switch(attr(x, "type"),
           "right"={
               temp <- x[,2]
               temp <- ifelse(is.na(temp), "?", ifelse(temp==0, "+"," "))
               paste(format(x[,1]), temp, sep='')
           },
           "counting"= {
               temp <- x[,3]
               temp <- ifelse(is.na(temp), "?", ifelse(temp==0, "+",""))
               paste('(', format(x[,1]), ',', format(x[,2]), temp,
                     ']', sep='')
           },
           "left" ={
               temp <- x[,2]
               temp <- ifelse(is.na(temp), "?", ifelse(temp==0, "<"," "))
               paste(temp, format(x[,1]), sep='')
           },
           "interval"= {
               stat <- x[,3]
               temp <- c("+", "", "-", "]")[stat+1]
               temp2 <- ifelse(stat==3,
			 paste("[", format(x[,1]), ", ",format(x[,2]), sep=''),
			 format(x[,1]))
               ifelse(is.na(stat), "NA", paste(temp2, temp, sep=''))
           },
           "mright" = {  #multi-state
               temp <- x[,2]
               end <- c("+", paste(":", attr(x, "states"), sep='')) #endpoint
               temp <- ifelse(is.na(temp), "?", end[temp+1])
               paste(format(x[,1]), temp, sep='')
           },
           "mcounting"= {
               temp <- x[,3]
               end <- c("+", paste(":", attr(x, "states"), sep='')) #endpoint
               temp <- ifelse(is.na(temp), "?", end[temp+1])
               paste('(', format(x[,1]), ',', format(x[,2]), temp,
                     ']', sep='')
           })
}


"[.Surv" <- function(x, i, j, drop=FALSE) {
    # If only 1 subscript is given, the result will still be a Surv object,
    #   and the drop argument is ignored.
    # I would argue that x[3:4,,drop=FALSE] should return a matrix, since
    #  the user has implicitly specified that they want a matrix.
    #  However, [.dataframe calls [.Surv with the extra comma; it's
    #  behavior drives the choice of default.
    if (missing(j)) {
        type <- attr(x, 'type')
        states <- attr(x, 'states')
        ctemp <- class(x)
        class(x) <- 'matrix'
        x <- x[i,, drop=FALSE]
        class(x) <- ctemp
        attr(x, 'type') <- type
        if (!is.null(states)) attr(x, "states") <- states
        x
    }
    else { # return  a matrix or vector
	if (is.R()) class(x) <- 'matrix'
        else oldClass(x) <- NULL
	NextMethod("[")
    }
}

is.na.Surv <- function(x) {
    as.vector(rowSums(is.na(unclass(x))) >0)
    }

Math.Surv <- function(...)  stop("Invalid operation on a survival time")
Ops.Surv  <- function(...)  stop("Invalid operation on a survival time")
Summary.Surv<-function(...) stop("Invalid operation on a survival time")
is.Surv <- function(x) inherits(x, 'Surv')
as.matrix.Surv <- function(x, ...) {
    y <- unclass(x)
    attr(y, "type") <- NULL
    attr(y, "states") <- NULL
    y
    }
