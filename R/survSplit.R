survSplit<-function(data, cut, end, event, start="tstart", id, zero=0,
                     episode) {
    if (!is.numeric(cut) || any(!is.finite(cut)))
        stop("cut must be a vector of finite numbers")
    cut<-sort(cut)
    ntimes <- length(cut)
    n <- nrow(data)
    if (missing(end) || missing(event))
        stop("the end and event arguments are required")

    if (is.character(event) && length(event) ==1 &&
        event %in% names(data)) status <- data[[event]]
    else stop("'event' must be a variable name in the data set")

    if (is.character(end) && length(end) ==1 &&
        end %in% names(data)) time2 <- data[[end]]
    else stop("'end' must be a variable name in the data set")

    if (!(is.character(start) && length(start)==1))
        stop("'start' must be a variable name")
    else {
        if (start %in% names(data)) time1 <- data[[start]]
        else {
            time1 <- rep(0., length.out = n)
            start <- make.names(start)  #force valid syntax
        }
    }
    if (any(time1 >= time2))
        stop("start time must be < stop time")

    if (!missing(id)) {
        if (!(is.character(id) && length(id)==1))
            stop("'id' must be a variable name")
        else {
            if (!(id %in% names(data)))
                data[[make.names(id)]] <- 1:n
        }
    }

    if (is.logical(status)) censor <- FALSE
    else if (is.factor(status)) censor <- levels(status)[1]
    else censor <- max(status) -1

    storage.mode(time1) <- "double"
    storage.mode(time2) <- "double"
    index <- .Call("survsplit", time1, time2, as.double(cut))
    newdata <- data[index$row,]
    row.names(newdata) <- NULL    # erase R's manufactured row names
    newdata[[start]] <- index$start
    newdata[[end]]   <- index$end
    status <- newdata[[event]]
    status[index$censor] <- censor
    newdata[[event]] <- status

    if (!missing(episode)) {
        if (!is.character(episode)) stop("episode must be a character string")
        newdata[[make.names(episode)]] <- index$interval +1
    }
    newdata
}

