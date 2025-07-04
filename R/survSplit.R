survSplit <- function(formula, data, subset, na.action=na.pass,
                              cut, start="tstart", id, zero=0, episode,
                              end="tstop", event="event", added) {
    Call <- match.call()
    if (missing(formula) || is.data.frame(formula)) {
        # an old style call
        # match arguments and build a formula
        if (missing(data)) {
            if (!missing(formula)) {
                names(Call)[[2]] <- "data"
                # The line above is sneaky: it makes model.frame() work later
                data <- formula
            }
            else  stop("a data frame is required")
        }
        if (missing(end) || missing(event))
            stop("either a formula or the end and event arguments are required")

        if (!(is.character(event) && length(event) ==1 &&
              event %in% names(data)))
            stop("'event' must be a variable name in the data set")

        if (!(is.character(end) && length(end) ==1 &&
              end %in% names(data)))
            stop("'end' must be a variable name in the data set")
    
        if (!(is.character(start) && length(start)==1))
            stop("'start' must be a variable name")
        if (start %in% names(data)) temp <- paste(start, end, event, sep=',')
        else temp <- paste(end, event, sep=',')
        
        formula <- as.formula(paste("Surv(", temp, ")~ ."))
    }
    else if (missing(formula)) 
        stop("either a formula or the end and event arguments are required")
 
    # create a call to model.frame() that contains the formula (required)
    #  and any other of the relevant optional arguments
    # then evaluate it in the proper frame
    indx <- match(c("data", "weights", "subset"),
                  names(Call), nomatch=0) 
    temp <- Call[c(1L,indx)]  # only keep the arguments we wanted
    temp$formula <- formula
    temp$na.action <- na.action
    temp[[1L]] <- quote(stats::model.frame)  # change the function called
    mf <- eval.parent(temp)      

    Y <- model.response(mf)
    states <- attr(Y, "states")
    if (!is.Surv(Y)) stop ("the model must have a Surv object as the response")
    if (!(attr(Y, "type") %in% c("right", "mright", "counting", "mcounting")))
        stop(paste("not valid for", attr(Y, "type"), "censored survival data"))
    nY <- ncol(Y)
    ymiss <- is.na(Y)  # these pass through unchanged
    if (nY ==2) {
        if (any(Y[!ymiss,1] <= zero))
            stop("'zero' parameter must be less than any observed times")
        Y <- cbind(zero, Y)
    }
    temp <- (Y[!ymiss,1] >= Y[!ymiss,2])
    if (any(temp)) stop("start time must be < stop time")
        
    if (!is.numeric(cut) || any(!is.finite(cut)))
        stop("cut must be a vector of finite numbers")
    cut <- unique(sort(cut))
    ntimes <- length(cut)
    n <- nrow(data)

    if (!missing(id)) {
        if (!is.character(id)) stop("id must be a variable name")
        if (id %in% names(mf)) stop("the suggested id name is already present")
        id <- make.names(id)
        if (id %in% names(mf)) stop("the suggested id name is already present")
        mf[[id]] <- 1:nrow(mf)
    }

    storage.mode(Y) <- "double"
    index <- .Call(Csurvsplit, Y[,1], Y[,2], as.double(cut))
    newdata <- mf[index$row, -1, drop=FALSE]
    row.names(newdata) <- NULL    # erase R's manufactured row names
    attr(newdata, "terms") <- NULL

    status <- Y[index$row, 3]
    status[index$censor] <- 0
    # the factor line needs 0:length(states) in case there is a state
    #  with no event, then status will be missing a category
    if (!is.null(states))  
        status <- factor(status, 0:length(states),labels=c("censor", states))

    # Did the user hand me a Surv call with multiple variables, or a
    #  premade Surv object?
    if (inherits(formula[[2]], "call") && formula[[2]][[1]]== as.name("Surv")){
        # it was a call, figure out the names
        # The user might have used something like Surv(status=abc, time=fred),
        #  so use match.call to find "abc" and "fred".  But give up if there
        #  is anything complex.
        temp <- match.call(Surv, formula[[2]])
        if (nY==2) {
            if (missing(end) && !is.null(temp[["time"]]) 
                && is.name(temp[["time"]]))
                end <- as.character(temp[["time"]])  # $time might match 'time2'
             
            if (missing(event) && !is.null(temp$time2) && is.name(temp$time2)) 
                event <- as.character(temp$time2)
            if (missing(event) && !is.null(temp$event) && is.name(temp$event))
                event <- as.character(temp$event)
        }
        else {
            if (missing(end) && !is.null(temp[["time"]]) 
                && is.name(temp["time"]))
                start <- as.character(temp[["time"]])
            if (missing(end) && !is.null(temp$time2) && is.name(temp$time2)) 
                end <- as.character(temp$time2)
            if (missing(event) && !is.null(temp$event) && is.name(temp$event))
                event <- as.character(temp$event)
            if (missing(start) && !is.null(temp$time) && is.name(temp$time))
                start <- as.character(temp$time)
        }

        newdata[[start]] <- index$start
        newdata[[end]]   <- index$end
        newdata[[event]] <- status
    }
    else {
        if (!inherits(formula[[2]], "name"))
            stop("left hand side not recognized")
        temp <- as.character(formula[[2]])
        newdata[temp] <- Surv(index$start, index$end, status)
    }

    if (!missing(episode)) {
        if (!is.character(episode)) stop("episode must be a character string")
        newdata[[make.names(episode)]] <- index$interval +1
    }
    if (!missing(added)) {
        if (!is.character(added)) stop("added must be a character string")
        newdata[[make.names(added)]] <- index$censor
    }
    newdata
}
